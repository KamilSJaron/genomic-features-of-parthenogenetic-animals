library('gsheet')
source('scripts/R_functions/load_genome_table.R')

tab_file <- "tables/assembly_table.tsv"
genome_tab <- load_genome_table()
genome_tab$haplotypes_assembled <- round(genome_tab[, 'assembly_size[M]'] / genome_tab[, 'haploid_length[M]'], 2)
genome_tab <- genome_tab[,c(1:2,5,6:12,19)]

# kick out all wo assembly
genome_tab <- genome_tab[!is.na(genome_tab[,'assembly_size[M]']),]

# read table with genome annotation stats
gene_annotations <- read.table('tables/gene_annotations.tsv', header = T, stringsAsFactors = F)

# add entries without annotation
to_add <- genome_tab$code[!genome_tab$code %in% gene_annotations$code]
indices_to_add <- nrow(gene_annotations) + 1:length(to_add)
gene_annotations[indices_to_add,] <- NA
gene_annotations[indices_to_add,'code'] <- to_add

# fix weird annotation format
row.names(gene_annotations) <- gene_annotations$code
gene_annotations['Dcor1','genes'] <- gene_annotations['Dcor1','mRNA']

genome_tab <- merge(gene_annotations[,c(1,2)], genome_tab)
row.names(genome_tab) <- genome_tab$code
# just sorting the columns
genome_tab <- genome_tab[, c(1, 3:12, 2)]

# genome_tab$gene_per_ploidy <- round(genome_tab$genes / genome_tab$ploidy)
genome_tab$gene_per_assembled_ploidy <- round(genome_tab$genes / genome_tab$haplotypes_assembled)

# add TEs
literature_data <- read.csv(text = gsheet2text("https://docs.google.com/spreadsheets/d/1T4BHQxzMGMlWiNJ7G9OLPXzdNBCqMw0mSO7C_ggGEbc/edit?usp=sharing", format='csv'),
                            stringsAsFactors = F, skip = 1, header = T, check.names = F)

extracted_TEs <- literature_data[,c('code','TEs [ % ]')]
colnames(extracted_TEs) <- c('code', 'annotated_TEs')

genome_tab <- merge(genome_tab, extracted_TEs)
genome_tab$annotated_TEs[is.na(genome_tab$annotated_TEs)] <- "NA"
rownames(genome_tab) <- genome_tab$code

# barplot(genome_tab$genes / genome_tab$haplotypes_assembled, ylim = c(0, 65000), col = 'green')
# barplot(genome_tab$genes / genome_tab$ploidy, ylim = c(0, 65000), add = T)

genome_tab <- genome_tab[, c(2,4:14)]
write.table(genome_tab, tab_file, quote = F, sep = '\t', row.names = F)


