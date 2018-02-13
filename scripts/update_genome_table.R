tab_file <- 'tables/genome_table.tsv'
genome_tab <- read.table(tab_file, header = T, stringsAsFactors = F)

stat_files <- list.files(path = ".", pattern = "genome.stats", recursive = T, include.dirs = T)
#	assembly_size number_of_scaffolds N50

check_sp_and_get_row <- function(sp){
	row <- sp == genome_tab$code
	if ( ! any(row) ){
		genome_tab[nrow(genome_tab)+1,] <- NA
		genome_tab[nrow(genome_tab), 'code'] <- sp
		row <- sp == genome_tab$code
	}
	return(row)
}

get_value <- function(line){
	as.numeric(strsplit(stat_lines[line], '\t')[[1]][[2]])
}

for(stat_file in stat_files){
	sp <- strsplit(stat_file, "/")[[1]][2]
	row <- check_sp_and_get_row(sp)

	stat_lines <- readLines(stat_file)

	genome_tab$assembly_size[row] <- get_value(1)
	genome_tab$number_of_scaffolds[row] <- get_value(2)
	genome_tab$N50[row] <- get_value(7)
}

TE_files <- list.files(path = "data", pattern = "Counts.txt", recursive = T, include.dirs = T)

for(TE_file in TE_files){
	sp <- strsplit(TE_file, "/")[[1]][1]
	row <- check_sp_and_get_row(sp)

	TEs <- read.table(paste0('data/',TE_file))

	genome_tab[row,'TEs'] <- round(sum(TEs[1:6,'V2']) / TEs[nrow(TEs),'V2'], 3) * 100
	genome_tab[row,'other_repeats'] <- round(sum(TEs[7:12,'V2']) / TEs[nrow(TEs),'V2'], 3) * 100
	genome_tab[row,'all_repeats'] <- round(sum(TEs[-nrow(TEs),'V2']) / TEs[nrow(TEs),'V2'], 3) * 100
}

write.table(genome_tab, tab_file, quote = F, sep = '\t', row.names = F)