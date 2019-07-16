#################
### load data ###
#################

library('gsheet')

# table to update
tab_file <- 'tables/genome_table.tsv'
genome_tab <- read.table(tab_file, header = T, stringsAsFactors = F, skip = 1, check.names = F)
row.names(genome_tab) <- genome_tab$code

# following line corrects names that are automatically replaced when R loads the table
# colnames(genome_tab)[c(3:5, 13)] <- c('assembly_size[M]', 'number_of_scaffolds[k]', 'N50[k]', 'haploid_length[M]')
# download table for cases if there is a new organism
dl_table <- read.table('tables/download_table.tsv', header = T, row.names = 1, stringsAsFactors = F)

sp_with_reads <- rownames(dl_table)[!is.na(dl_table$reads)]
sp_with_genomes <- rownames(dl_table)[!is.na(dl_table$genome)]

#####################################################
#                   COMPUTED DATA                   #
# ASSEMBLY STATS ( scripts/fasta2genomic_stats.py ) #
#    assembly_size number_of_scaffolds N50          #
#####################################################

# funciton to check files
source('scripts/R_functions/checkFiles.R')
# functions to raed GenomeScope and BUSCO outputs
source('scripts/R_functions/output_parsers.R')

expand_table_if_needed <- function(.sp, .genome_tab){
    row <- .sp == .genome_tab$code
    row[is.na(row)] <- FALSE
    if ( ! any(row) ){
        row <- nrow(.genome_tab) + 1
        .genome_tab[row, ] <- NA
        .genome_tab[row, 'code'] <- .sp
        # fill full latin names from download table : 'tables/download_table.tsv'
        .genome_tab[row, 'species'] <- dl_table[sp,'species']
    }
    return(.genome_tab)
}

stat_files <- paste("data", sp_with_genomes, "genome.stats", sep = "/")
stat_files <- checkFiles(stat_files, 'stats')

for(stat_file in stat_files){
    sp <- ssplit(stat_file, "/")[2]
    genome_tab <- expand_table_if_needed(sp, genome_tab)
    row <- sp == genome_tab$code

    stat_lines <- readLines(stat_file)

    genome_tab[row, 'assembly_size[M]'] <- round(get_value(1) / 1000000, 1)
    genome_tab[row, 'number_of_scaffolds[k]'] <- round(get_value(2) / 1000, 1)
    genome_tab[row, 'N50[k]'] <- round(get_value(7) / 1000, 1)
}

#############################################
# dnaPipeTE ( scripts/annotate_repeats.sh ) #
# TEs other_repeats all_repeats             #
#############################################

TE_files <- paste("data", sp_with_reads, "dnaPipeTE/Counts.txt", sep = "/")
TE_files <- checkFiles(TE_files, 'TEs')

for(TE_file in TE_files){
    sp <- ssplit(TE_file, "/")[2]
    genome_tab <- expand_table_if_needed(sp, genome_tab)
    row <- sp == genome_tab$code

    TEs <- read.table(TE_file)

    genome_tab[row,'TEs'] <- round(sum(TEs[1:6,'V2']) / TEs[nrow(TEs),'V2'], 3) * 100
    genome_tab[row,'other_repeats'] <- round(sum(TEs[7:12,'V2']) / TEs[nrow(TEs),'V2'], 3) * 100
    genome_tab[row,'all_repeats'] <- round(sum(TEs[-nrow(TEs),'V2']) / TEs[nrow(TEs),'V2'], 3) * 100
}

#############################################
# BUSCO ( scripts/busco.sh )                #
# complete duplicated fragmented missing    #
#############################################

busco_files <- paste("data", sp_with_genomes, "busco/short_summary_busco.txt", sep = "/")
busco_files <- checkFiles(busco_files, 'busco')

for(busco_file in busco_files){
    sp <- strsplit(busco_file, "/")[[1]][2]
    genome_tab <- expand_table_if_needed(sp, genome_tab)
    row <- sp == genome_tab$code

    genome_tab[row, c('complete', 'fragmented', 'duplicated', 'missing')] <- read_busco(busco_file)
}

###################################################
# Smudgeplot ( scripts/generate_smudgeplot.sh )   #
#            Ploidy                               #
###################################################

smudgeplot_files <- paste0("data/", sp_with_reads, "/smudgeplot/", sp_with_reads, "_verbose_summary.txt")
smudgeplot_files <- checkFiles(smudgeplot_files, 'smudgeplot files')

for(smudge_file in smudgeplot_files){
    sp <- ssplit(smudge_file, "/")[2]
    genome_tab <- expand_table_if_needed(sp, genome_tab)
    row <- sp == genome_tab$code

    smudge_est_ploidy <- ssplit(readLines(smudge_file)[6], '\t')[2]

    genome_tab[row,'ploidy'] <- as.numeric(smudge_est_ploidy)
}

# This species is most likely hozmozygous within genome copies,
# but very heterozygous between othnologs, it is AABB where AA and BB are the same,
# but AB is very diverged (~12%)
genome_tab[genome_tab$code %in% c("Avag1", "Rmag1"),'ploidy'] <- 4

#############################################################
# GenomeScope ( scripts/GenomeScope.sh )                    #
# kmer_genome_size, kmer_heterozygosity, kmer_repetitions   #
#############################################################

# genome_tab[, c('haploid_length[M]', 'repeats', 'heterozygosity')] <- NA
genomescope_files <- paste("data", sp_with_reads, "genomescope/summary.txt", sep = "/")
genomescope_files <- checkFiles(genomescope_files, 'genomescope files (not converged)')

for(genomescope_file in genomescope_files){
    sp <- ssplit(genomescope_file, "/")[2]
    genome_tab <- expand_table_if_needed(sp, genome_tab)
    row <- sp == genome_tab$code

    genome_tab[row, c('haploid_length[M]',
                      'repeats',
                      'heterozygosity')] <- parse_genomescope_summary(genomescope_file)
}

# tetraploids have for some reason a problem with summaries, all other ploidies have values corresponding to figures, tetraplods do not
# therefore I use values from their figures:

genome_tab[c('Aric1', 'Rmac1', 'Rmag1', 'Mjav2', 'Mare2'),
           'heterozygosity'] <- c(6.1, 12.4, 11.7, 8.5, 8.1)

#####################################################
###                LITERATURE DATA                ###
### add reproduction mode, ploidy and genome_size ###
#####################################################

literature_data <- read.csv(text = gsheet2text("https://docs.google.com/spreadsheets/d/1T4BHQxzMGMlWiNJ7G9OLPXzdNBCqMw0mSO7C_ggGEbc/edit?usp=sharing", format='csv'),
                            stringsAsFactors = F, skip = 1, header = T, check.names = F)
literature_data <- literature_data[,c(2, 8, 9)]

repr_col <- 'cellular mechanism of asexuality'

mitotic <- grepl('mitotic', literature_data[,repr_col]) | grepl('endoduplication', literature_data[,repr_col])
central_fusion <- grepl('central', literature_data[,repr_col])
terminal_fusion <- grepl('terminal', literature_data[,repr_col])
gdupl <- grepl('gamete dupl', literature_data[,repr_col])
unknown_meiotic <- grepl('meiotic', literature_data[,repr_col]) & ! (central_fusion | terminal_fusion | mitotic)

# create three categories
literature_data$callular_mechanism <- NA # unknown is default
literature_data$callular_mechanism[mitotic] <- 'functional_apomixis' # mitosis or endoduplication
literature_data$callular_mechanism[central_fusion] <- 'central_fusion'  # automixis central fusion
literature_data$callular_mechanism[terminal_fusion] <- 'terminal_fusion' # automixis terminal fusion
literature_data$callular_mechanism[gdupl] <- 'gamete_duplication'  # gamete duplications
literature_data$callular_mechanism[unknown_meiotic] <- 'unknown_automixis' # unknown automixis

literature_data <- literature_data[literature_data$code %in% genome_tab$code,]

genome_tab$callular_mechanism <- NA
genome_tab$hybrid_origin <- NA

row.names(genome_tab) <- genome_tab$code
columns <- c('callular_mechanism', 'hybrid_origin')
genome_tab[literature_data$code, columns] <- literature_data[, columns]

######################
### Sort the table ###
######################

# soring rows
# lines
#  vertebrates
#  rotifers
#  arthropods
#  nematodes
#  tardigrades
desired_order <- c('Pfor1',
                   'Avag1', 'Aric1', 'Rmac1', 'Rmag1',
                   'Lcla1', 'Tpre1', 'Obir1', 'Amel1', 'Amel2', 'Amel3', 'Aruf1', 'Fcan1', 'Dpul1', 'Dpul2', 'Dpul3', 'Pvir1',
                   'Psam1', 'Mbel1', 'Dcor1', 'Dpac1', 'Pdav1', 'Ps591', 'Ps791', 'Minc1', 'Minc2', 'Mjav1', 'Mjav2', 'Mare1', 'Mare2', 'Mare3', 'Mflo1', 'Ment1', 'Anan1',
                   'Hduj1', 'Rvar1')
if ( length(desired_order) == nrow(genome_tab) ){
    genome_tab <- genome_tab[desired_order, ]
} else {
    cat('There are species in the table with unknown position by the hand-made-order in this script.\n')
    cat('namely : ')
    cat( genome_tab$code[ !genome_tab$code %in% desired_order] )
}

# sorting columns
genome_tab <- genome_tab[, c('code', 'species', 'callular_mechanism', 'hybrid_origin', 'ploidy',
                             'assembly_size[M]', 'number_of_scaffolds[k]', 'N50[k]',
                             'complete', 'fragmented', 'duplicated', 'missing',
                             'haploid_length[M]', 'heterozygosity', 'repeats',
                             'TEs','other_repeats','all_repeats')]

######################

extra_header <- c(rep('-', 2),
                  'literature', rep('-', 1),
                  'smudgeplot',
                  'assembly', rep('-', 2),
                  'BUSCO', rep('-', 3),
                  'GenomeScope', rep('-', 2),
                  'dnaPipeTE', rep('-', 2))
asm_template <- matrix(ncol = length(extra_header))
colnames(asm_template) <- extra_header
header <- as.data.frame(asm_template)[FALSE, ]

write.table(header, tab_file, quote = F, sep = '\t', row.names = F)
write.table(genome_tab, tab_file, quote = F, sep = '\t', row.names = F, append = T)
