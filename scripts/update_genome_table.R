#################
### load data ###
#################

library('gsheet')

# table to update
tab_file <- 'tables/genome_table.tsv'
genome_tab <- read.table(tab_file, header = T, stringsAsFactors = F, skip = 1, check.names = F)
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

checkFiles <- function(.files, .type){
    files_exists <- file.exists(.files)
    if( !all(files_exists) ){
        warning(paste("missing ", .type, "files:\n", paste('\t\t\t\t', .files[!files_exists], collapse = '\n')))
        .files <- .files[files_exists]
    }
    return(.files)
}

stat_files <- paste("data", sp_with_genomes, "genome.stats", sep = "/")
stat_files <- checkFiles(stat_files, 'stats')

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

get_value <- function(line){
    as.numeric(strsplit(stat_lines[line], '\t')[[1]][[2]])
}

# function that unlists strsplit
ssplit <- function (s, split = "="){
    unlist(strsplit(s, split = split))
}

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

read_busco <- function(busco_file){
    if( is.na(busco_file) ){
        return( rep(NA, 4) )
    }
    busco_file <- readLines(busco_file)
    total_genes <- as.numeric(ssplit(busco_file[15], '\t')[2])
    bscores <- c(complete = as.numeric(ssplit(busco_file[10], '\t')[2]),
                 fragmented = as.numeric(ssplit(busco_file[13], '\t')[2]),
                 duplicated = as.numeric(ssplit(busco_file[12], '\t')[2]),
                 missing = as.numeric(ssplit(busco_file[14], '\t')[2]))
    bscores <- round(100 * (bscores / total_genes), 2)
    return(bscores)
}

busco_files <- paste("data", sp_with_genomes, "busco/short_summary_busco.txt", sep = "/")
busco_files <- checkFiles(busco_files, 'busco')

for(busco_file in busco_files){
    sp <- strsplit(busco_file, "/")[[1]][2]
    genome_tab <- expand_table_if_needed(sp, genome_tab)
    row <- sp == genome_tab$code

    genome_tab[row, c('complete', 'fragmented', 'duplicated', 'missing')] <- read_busco(busco_file)
}

####################################
# ATLAS ( scripts/est_theta.sh )   #
# ML_heterozygosity                #
####################################

# TODO

#############################################################
# GenomeScope ( scripts/GenomeScope.sh )                    #
# kmer_genome_size, kmer_heterozygosity, kmer_repetitions   #
#############################################################

parse_genomescope_summary <- function(file){
    genoscope_file <- readLines(file)
    line <- ssplit(genoscope_file[5], ' ')
    het_string <- line[line != ''][3]
    heterozygosity <- as.numeric(substr(het_string, 0, nchar(het_string) - 1))

    line <- ssplit(genoscope_file[6], ' ')
    haploid_genome <- as.numeric(gsub(",", "", line[line != ''][6]))

    line <- ssplit(genoscope_file[7], ' ')
    repeats <- (as.numeric(gsub(",", "", line[line != ''][6])) * 100) / haploid_genome

    if( heterozygosity == -100 ){
        return( c(NA, NA, NA) )
    } else {
        return( c(round(haploid_genome / 1e6, 1), round(repeats, 2), round(heterozygosity, 2)) )
    }
}

polyploid_genomescope_files <- paste("data", sp_with_reads, "genomescope_v2/summary.txt", sep = "/")
polyploid_genomescope_files <- checkFiles(polyploid_genomescope_files, 'polyploid genomescope v2')

genomescope_files <- paste("data", sp_with_reads[!sp_with_reads %in% substring(polyploid_genomescope_files, 6, 10)], "genomescope/summary.txt", sep = "/")
genomescope_files <- checkFiles(genomescope_files, 'genomescope')

genomescope_files <- checkFiles(genomescope_files, 'genomescope')
genomescope_files <- c(polyploid_genomescope_files, genomescope_files)

for(genomescope_file in genomescope_files){
    sp <- ssplit(genomescope_file, "/")[2]
    genome_tab <- expand_table_if_needed(sp, genome_tab)
    row <- sp == genome_tab$code

    genome_tab[row, c('haploid_length[M]', 'repeats', 'heterozygosity')] <- parse_genomescope_summary(genomescope_file)
}

#####################################################
###                LITERATURE DATA                ###
### add reproduction mode, ploidy and genome_size ###
#####################################################

literature_data <- read.csv(text = gsheet2text("https://docs.google.com/spreadsheets/d/1T4BHQxzMGMlWiNJ7G9OLPXzdNBCqMw0mSO7C_ggGEbc/edit?usp=sharing", format='csv'),
                            stringsAsFactors = F, skip = 1, header = T, check.names = F)
literature_data <- literature_data[,c(2, 8, 9, 12)]

meio <- grepl('auto', literature_data[,'reproduction_mode']) | grepl('meiosis', literature_data[,'reproduction_mode'])
mito <- grepl('apo', literature_data[,'reproduction_mode']) | grepl('mitosis', literature_data[,'reproduction_mode'])
gupl <- grepl('dupl', literature_data[,'reproduction_mode'])

# create three categories
literature_data$reproduction_mode[meio] <- 'automixis' # involve meiosis
literature_data$reproduction_mode[mito] <- 'apomixis' # wo mitosis
literature_data$reproduction_mode[gupl] <- 'gamete_duplication' # gamete duplications
literature_data$reproduction_mode[ ! (meio | mito | gupl) ] <- NA # unknown

literature_data <- literature_data[literature_data$code %in% genome_tab$code,]

genome_tab$ploidy <- NA
genome_tab$reproduction_mode <- NA
genome_tab$hybrid_origin <- NA

row.names(genome_tab) <- genome_tab$code
columns <- c('ploidy', 'reproduction_mode', 'hybrid_origin')
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
                   'Lcla1', 'Tpre1', 'Obir1', 'Aruf1', 'Fcan1', 'Dpul1', 'Dpul2', 'Dpul3', 'Dpul4', 'Dpul5', 'Pvir1',
                   'Psam1', 'Dcor1', 'Dpac1', 'Pdav1', 'Ps591', 'Ps791', 'Minc1', 'Minc2', 'Minc3', 'Mjav1', 'Mjav2', 'Mare1', 'Mare2', 'Mare3', 'Mflo1', 'Mflo2', 'Ment1', 'Anan1',
                   'Hduj1', 'Rvar1')
if ( length(desired_order) == nrow(genome_tab) ){
    genome_tab <- genome_tab[desired_order, ]
} else {
    cat('There are species in the table with unknown position by the hand-made-order in this script.\n')
    cat('namely : ')
    cat( genome_tab$code[ !genome_tab$code %in% desired_order] )
}

# sorting columns
genome_tab <- genome_tab[, c('code', 'species', 'reproduction_mode', 'hybrid_origin', 'ploidy',
                             'assembly_size[M]', 'number_of_scaffolds[k]', 'N50[k]',
                             'complete', 'fragmented', 'duplicated', 'missing',
                             'haploid_length[M]', 'heterozygosity', 'repeats',
                             'TEs','other_repeats','all_repeats')]

######################

extra_header <- c(rep('-', 5), 'assembly', rep('-', 2), 'BUSCO', rep('-', 3), 'GenomeScope', rep('-', 2), 'dnaPipeTE', rep('-', 2) )
asm_template <- matrix(ncol = length(extra_header))
colnames(asm_template) <- extra_header
header <- as.data.frame(asm_template)[FALSE, ]

write.table(header, tab_file, quote = F, sep = '\t', row.names = F)
write.table(genome_tab, tab_file, quote = F, sep = '\t', row.names = F, append = T)