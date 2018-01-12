tab_file <- 'tables/genome_table.tsv'
genome_tab <- read.table(tab_file, header = T)

stat_files <- list.files(path = ".", pattern = "genome.stats", recursive = T, include.dirs = T)
#	assembly_size number_of_scaffolds N50

for(stat_file in stat_files){
	sp <- strsplit(stat_file, "/")[[1]][2]
	stat_lines <- readLines(stat_file)

	genome_tab$assembly_size[sp == genome_tab$code] <- as.numeric(strsplit(stat_lines[1], '\t')[[1]][[2]])
	genome_tab$number_of_scaffolds[sp == genome_tab$code] <- as.numeric(strsplit(stat_lines[2], '\t')[[1]][[2]])
	genome_tab$N50[sp == genome_tab$code] <- as.numeric(strsplit(stat_lines[7], '\t')[[1]][[2]])
}

write.table(genome_tab, tab_file, quote = F, sep = '\t', row.names = F)