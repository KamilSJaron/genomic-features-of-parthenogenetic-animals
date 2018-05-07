tab_file <- 'tables/genome_table.tsv'
genome_tab <- read.table(tab_file, header = T, stringsAsFactors = F, skip = 1)
colnames(genome_tab)[3:5] <- c('assembly_size[M]', 'number_of_scaffolds[k]', 'N50[k]')
dl_table <- read.table('tables/download_table.tsv', header = T, row.names = 1, stringsAsFactors = F)

#####################################################
# ASSEMBLY STATS ( scripts/fasta2genomic_stats.py ) #
#	assembly_size number_of_scaffolds N50           #
#####################################################

stat_files <- list.files(path = ".", pattern = "genome.stats", recursive = T, include.dirs = T)

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

# TODO

####################################
# ATLAS ( scripts/est_theta.sh )   #
# ML_heterozygosity                #
####################################

#############################################################
# GenomeScope ( scripts/GenomeScope.sh )                    #
# kmer_genome_size, kmer_heterozygosity, kmer_repetitions   #
#############################################################

write.table(genome_tab, tab_file, quote = F, sep = '\t', row.names = F)