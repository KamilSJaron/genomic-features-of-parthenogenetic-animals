args <- commandArgs(trailingOnly=TRUE)
sp <- args[1]

# USING KMERS
genome_tab <- read.table('tables/genome_table.tsv', skip = 1, header = T, check.names = F)

my_sp <- genome_tab$code == sp
haploid_genome <- genome_tab[my_sp, 'haploid_length[M]'] * 1e6

if( all(is.na(haploid_genome)) ){
    # USING ASSEMBLY
    genome_table <- read.table('tables/genome_table.tsv', header = T, skip = 1)
    # extract full species name from table
    sp_name <- genome_table[grepl(substr(sp, 1, 4), genome_table$code), 'species']
    # identify all genomes of species
    asm_sizes <- genome_table[genome_table$species == sp_name,'assembly_size.M.'] * 1e6
    # print the size of longest
    cat(sprintf("%.0f",max(asm_sizes, na.rm = T)))
} else {
    cat(sprintf("%.0f", haploid_genome))
}
