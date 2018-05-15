args <- commandArgs(trailingOnly=TRUE)
sp <- args[1]

# load table
genome_table <- read.table('tables/genome_table.tsv', header = T, skip = 1)
# extract full species name from table
sp_name <- genome_table[grepl(sp, genome_table$code), 'species']
# identify all genomes of species
asm_sizes <- genome_table[genome_table$species == sp_name,'assembly_size.M.']
# print the size of longest
cat(max(asm_sizes, na.rm = T))
