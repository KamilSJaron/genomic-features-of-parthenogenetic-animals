source('scripts/R_functions/load_genome_table.R')

tab_file <- "tables/genome_table_infered_from_reads.tsv"

genome_tab <- load_genome_table()
asm_sizes <- genome_tab[c(22, 33),'assembly_size[M]']

genome_tab <- genome_tab[,c(2,5,13:17)]

samples_with_values <- !apply(is.na(genome_tab[,3:7]), 1, all)

genome_tab <- genome_tab[samples_with_values,]

names_split <- strsplit(genome_tab$species, "_")
genus_names <- sapply(genome_tab$species, substr, 1, 1)
species_names <- sapply(names_split, function(x){ x[2] })
sp_labels <- paste(genus_names, species_names, sep = '. ')
genome_tab$species <- sp_labels

write.table(genome_tab, tab_file, quote = F, sep = '\t', row.names = F)