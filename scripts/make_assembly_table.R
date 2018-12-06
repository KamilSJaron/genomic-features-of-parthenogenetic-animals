source('scripts/R_functions/load_genome_table.R')

tab_file <- "tables/assembly_table.tsv"
genome_tab <- load_genome_table()
genome_tab$haplotypes_assembled <- round(genome_tab[, 'assembly_size[M]'] / genome_tab[, 'haploid_length[M]'], 2)
genome_tab <- genome_tab[,c(1:2,5,6:12,19)]

gene_annotations <- read.table('tables/gene_annotations.tsv', header = T)
gene_annotations[c(28, 29),] <- NA
gene_annotations$code <- as.character(gene_annotations$code)
gene_annotations[c(28, 29),'code'] <- c("Mare3", "Aruf1")
row.names(gene_annotations) <- gene_annotations$code
gene_annotations['Dcor1','genes'] <- gene_annotations['Dcor1','mRNA']

genome_tab <- merge(genome_tab, gene_annotations[,c(1,2)])

row.names(genome_tab) <- genome_tab$code
sp_codes <- load_genome_table(1:4)[,1]
sp_codes <- sp_codes[sp_codes %in% genome_tab$code]

# genome_tab$gene_per_ploidy <- round(genome_tab$genes / genome_tab$ploidy)
genome_tab$gene_per_assembled_ploidy <- round(genome_tab$genes / genome_tab$haplotypes_assembled)

# barplot(genome_tab$genes / genome_tab$haplotypes_assembled, ylim = c(0, 65000), col = 'green')
# barplot(genome_tab$genes / genome_tab$ploidy, ylim = c(0, 65000), add = T)

genome_tab <- genome_tab[sp_codes, c(2,4:14)]
write.table(genome_tab, tab_file, quote = F, sep = '\t', row.names = F)


