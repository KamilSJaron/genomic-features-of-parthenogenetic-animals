source('scripts/R_functions/checkFiles.R')

tab_file <- 'tables/palindrome_table.tsv'
species <- dir('data/', pattern = "*[0-9]")
filenames <- paste0('data/', species, '/MCScanX_biggap/', species, '_prot_palindrome_summary.txt')
filenames <- checkFiles(filenames, "palindrome summary file")

species <- substr(filenames, 6, 10)

palindorme_files <- paste0('data/', species, '/MCScanX_biggap/', species, '_prot_palindromes.collinearity')
palindorme_files <- checkFiles(palindorme_files, "palindrome collinearity file")

process_palindorme_file <- function(.x){
    genome_tab <- read.table(.x, skip = 5, header = T)
    return(c(sum(genome_tab$spacer < 170000), nrow(genome_tab), round(median(genome_tab$spacer)), range(genome_tab$spacer)))
}

which_are_close <- function(.x){
    genome_tab <- read.table(.x, skip = 5, header = T)
    return(genome_tab$spacer < 170000)
}

fish_out_reverse_blocks <- function(.x){
    as.numeric(strsplit(readLines(.x)[3], '\t ')[[1]][2])
}

get_palindrome_sizes <- function(.x){
    lines <- readLines(.x)
    aligments <- rep(0, length(lines))
    aligments[grepl("Alig", lines)] <- 1
    2 * (table(cumsum(aligments)) - 1)
}

palindorme_summaries <- lapply(filenames, process_palindorme_file)
reverse_blocks <- sapply(filenames, fish_out_reverse_blocks)

keep <- reverse_blocks > 0
palindromes <- sapply(palindorme_summaries[keep], function(x){ x[2] })
close_palindromes <- sapply(palindorme_summaries[keep], function(x){ x[1] })
palindrome_genes <- lapply(palindorme_files[keep], get_palindrome_sizes)
which_close_palindromes <- lapply(filenames[keep], which_are_close)

potentially_affected_genes <- sapply(palindrome_genes, sum)
affected_genes <- sapply(1:length(palindrome_genes), function(x){ sum(palindrome_genes[[x]][which_close_palindromes[[x]]]) })

pal_tab <- data.frame(code = species[keep], reverse_blocks = reverse_blocks[keep], palindromes = palindromes, close_palindromes = close_palindromes, potentially_affected_genes = potentially_affected_genes, affected_genes = affected_genes)

gene_annotations <- read.table('tables/gene_annotations.tsv', header = T)
source('scripts/R_functions/load_genome_table.R')
species_names <- load_genome_table(1:4)[,c(1,2)]

pal_tab <- merge(merge(pal_tab, gene_annotations), species_names)
names_split <- strsplit(pal_tab$species, "_")
genus_names <- sapply(names_split, substr, 1, 1)[1,]
species_names <- sapply(names_split, function(x){ x[2] })
sp_labels <- paste(genus_names, species_names, sep = '. ')

pal_tab <- data.frame(species = sp_labels, reverse_blocks = reverse_blocks[keep], palindromes = palindromes, cloese_palindromes = close_palindromes, potentially_affected_genes = potentially_affected_genes, affected_genes = affected_genes, annotated_genes = pal_tab$genes, propotion_of_affected_genes = round((pal_tab$potentially_affected_genes / pal_tab$genes) * 100, 2))
row.names(pal_tab) <- 1:9

pal_tab <- pal_tab[c(8, 1, 7, 3, 2, 6, 5, 4, 9),]

write.table(pal_tab, tab_file, quote = F, sep = '\t', row.names = F)
