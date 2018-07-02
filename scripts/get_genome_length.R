args <- commandArgs(trailingOnly=TRUE)
sp <- args[1]

# USING KMERS
genomescope_file <- paste('data', sp, 'genomescope/summary.txt', sep = '/')

genomescope_genome_size_line <- unlist(strsplit(readLines(genomescope_file)[6], ' '))
haploid_genome <- as.numeric(gsub(",", "", genomescope_genome_size_line[genomescope_genome_size_line != ''][6]))

if( haploid_genome == -1 | is.na(haploid_genome) ){
    # USING ASSEMBLY
    genome_table <- read.table('tables/genome_table.tsv', header = T, skip = 1)
    # extract full species name from table
    sp_name <- genome_table[grepl(substr(sp, 1, 4), genome_table$code), 'species']
    # identify all genomes of species
    asm_sizes <- genome_table[genome_table$species == sp_name,'assembly_size.M.'] * 1e6
    # print the size of longest
    cat(max(asm_sizes, na.rm = T))
} else {
    cat(haploid_genome)
}