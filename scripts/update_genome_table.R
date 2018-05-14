#################
### load data ###
#################

# table to update
tab_file <- 'tables/genome_table.tsv'
genome_tab <- read.table(tab_file, header = T, stringsAsFactors = F, skip = 1)
# following line corrects names that are automatically replaced when R loads the table
colnames(genome_tab)[c(3:5, 13)] <- c('assembly_size[M]', 'number_of_scaffolds[k]', 'N50[k]', 'haploid_length[M]')
# download table for cases if there is a new organism
dl_table <- read.table('tables/download_table.tsv', header = T, row.names = 1, stringsAsFactors = F)

#####################################################
# ASSEMBLY STATS ( scripts/fasta2genomic_stats.py ) #
#	assembly_size number_of_scaffolds N50           #
#####################################################

stat_files <- list.files(path = "data", pattern = "genome.stats", recursive = T, include.dirs = T)
stat_files <- paste('data', stat_files, sep = '/')

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

for(stat_file in stat_files){
	sp <- strsplit(stat_file, "/")[[1]][2]
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

TE_files <- list.files(path = "data", pattern = "Counts.txt", recursive = T, include.dirs = T)

for(TE_file in TE_files){
	sp <- strsplit(TE_file, "/")[[1]][1]
	genome_tab <- expand_table_if_needed(sp, genome_tab)
	row <- sp == genome_tab$code

	TEs <- read.table(paste0('data/',TE_file))

	genome_tab[row,'TEs'] <- round(sum(TEs[1:6,'V2']) / TEs[nrow(TEs),'V2'], 3) * 100
	genome_tab[row,'other_repeats'] <- round(sum(TEs[7:12,'V2']) / TEs[nrow(TEs),'V2'], 3) * 100
	genome_tab[row,'all_repeats'] <- round(sum(TEs[-nrow(TEs),'V2']) / TEs[nrow(TEs),'V2'], 3) * 100
}

#############################################
# BUSCO ( scripts/busco.sh )                #
# complete duplicated fragmented missing    #
#############################################

ssplit <- function (s, split = "="){
    unlist(strsplit(s, split = split))
}

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

busco_files <- list.files(path = ".", pattern = "short_summary", recursive = T, include.dirs = T)

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

genomescope_files <- list.files(path = "data", pattern = "summary.txt", recursive = T, include.dirs = T)
genomescope_files <- paste('data', genomescope_files, sep = '/')

for(genomescope_file in genomescope_files){
	sp <- strsplit(genomescope_file, "/")[[1]][2]
	genome_tab <- expand_table_if_needed(sp, genome_tab)
	row <- sp == genome_tab$code

	genome_tab[row, c('haploid_length[M]', 'repeats', 'heterozygosity')] <- parse_genomescope_summary(genomescope_file)
}

######################
### Sort the table ###
######################

# soring rows
desired_order <- c('Pfor1',
                   'Avag1', 'Aric1', 'Rmac1', 'Rmag1',
				   'Lcla1', 'Cbir1', 'Fcan1',
				   'Pvir1', 'Dpul1',
				   'Pdav1', 'Dcor1', 'Dpac1', 'Minc1', 'Minc2', 'Minc3', 'Mjav1', 'Mjav2', 'Mare1', 'Mare2', 'Mflo1', 'Mflo2', 'Ment1',
				   'Hduj1', 'Rvar1')
if ( length(desired_order) == nrow(genome_tab) ){
	row.names(genome_tab) <- genome_tab$code
	genome_tab <- genome_tab[desired_order, ]
} else {
	cat('There are species in the table with unknown position by the hand-made-order in this script.\n')
	cat('namely : ')
	cat( genome_tab$code[ !genome_tab$code %in% desired_order] )
}

# sorting columns
genome_tab <- genome_tab[, c('code', 'species',
                             'assembly_size[M]', 'number_of_scaffolds[k]', 'N50[k]',
                             'complete', 'fragmented', 'duplicated', 'missing',
                             'haploid_length[M]', 'heterozygosity', 'repeats',
                             'TEs','other_repeats','all_repeats')]

######################

extra_header <- c(rep('-', 2), 'assembly', rep('-', 2), 'BUSCO', rep('-', 3), 'GenomeScope', rep('-', 2), 'dnaPipeTE', rep('-', 2) )
asm_template <- matrix(ncol = length(extra_header))
colnames(asm_template) <- extra_header
header <- as.data.frame(asm_template)[FALSE, ]

write.table(header, tab_file, quote = F, sep = '\t', row.names = F)
write.table(genome_tab, tab_file, quote = F, sep = '\t', row.names = F, append = T)