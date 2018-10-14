suppressWarnings(library('gsheet'))

literature_data <- read.csv(text = gsheet2text("https://docs.google.com/spreadsheets/d/1T4BHQxzMGMlWiNJ7G9OLPXzdNBCqMw0mSO7C_ggGEbc/edit?usp=sharing", format='csv'),
                            stringsAsFactors = F, skip = 1, header = T, check.names = F)


tab_file <- 'tables/genome_table.tsv'
genome_tab <- read.table(tab_file, header = T, stringsAsFactors = F, skip = 1, check.names = F)

# get unique species
genome_tab <- genome_tab[!genome_tab$code %in% c('Ps591', 'Ps791', 'Dpul1', 'Dpul2', 'Dpul3', 'Dpul4', 'Mjav1', 'Mare1', 'Mare3', 'Minc1'),]

diploids <- genome_tab$ploidy == 2 & !is.na(genome_tab$ploidy)
otherploids <- genome_tab$ploidy != 2 & !is.na(genome_tab$ploidy)

number_of_species <- length(unique(genome_tab$species))

# # genomes
cat("Number of reanalyzed genomes: ", number_of_species, "\n")

# genome sizes
cat("Range of haploid sizes [Mbp]:", paste(range(c(genome_tab[,'assembly_size[M]'], genome_tab[,'haploid_length[M]']), na.rm=T)), "\n")
cat("Haploid number of chromosomes: ", paste(range(literature_data$chromosomes_haploid, na.rm = T)), "\n")


# ploidy
cat("Number of diploid species: ", length(unique(genome_tab$species[diploids])), "\n")
cat("Number of polyploid species: ", length(unique(genome_tab$species[otherploids])), "\n")

cat("Max missing conserved genes in Meloidogyne: ", max(genome_tab[substr(genome_tab$species, 1, 11) == "Meloidogyne",'missing']), "\n")
cat("HGT in A. vaga: ", literature_data[literature_data$code == "Avag1", "HGT"], "\n")
cat("Heterozygosity of D. pachys: ", genome_tab[genome_tab$code == "Dpac1", 'heterozygosity'], "\n")

