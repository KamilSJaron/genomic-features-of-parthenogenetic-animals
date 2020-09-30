suppressWarnings(library('gsheet'))

literature_data <- read.csv(text = gsheet2text("https://docs.google.com/spreadsheets/d/1T4BHQxzMGMlWiNJ7G9OLPXzdNBCqMw0mSO7C_ggGEbc/edit?usp=sharing", format='csv'),
                            stringsAsFactors = F, skip = 1, header = T, check.names = F)


tab_file <- 'tables/genome_table.tsv'
genome_tab <- read.table(tab_file, header = T, stringsAsFactors = F, skip = 1, check.names = F)

# get unique species
genome_tab <- genome_tab[!genome_tab$code %in% c('Ps591', 'Ps791', 'Dpul1', 'Dpul4', 'Mjav1', 'Mare1', 'Mare3', 'Minc1'),]

diploids <- genome_tab$ploidy == 2 & !is.na(genome_tab$ploidy)
otherploids <- genome_tab$ploidy != 2 & !is.na(genome_tab$ploidy)

# genome_tab$heterozygosity[genome_tab$code %in% c('Avag1', 'Aric1', 'Rmac1', 'Rmag1')] <- 33.21
# range(genome_tab$heterozygosity[diploids], na.rm = T)
# range(genome_tab$heterozygosity[otherploids], na.rm = T)
# median(genome_tab$heterozygosity[diploids], na.rm = T)
# median(genome_tab$heterozygosity[otherploids], na.rm = T)

number_of_species <- length(unique(genome_tab$species))

genome_tab$callular_mechanism[is.na(genome_tab$callular_mechanism)] <- "unknown"
genome_tab$hybrid_origin[is.na(genome_tab$hybrid_origin)] <- "unknown"
cellular_mechs <- aggregate(callular_mechanism ~ species, data=genome_tab, FUN=function(x){x[1]})

mitotic_parth <- genome_tab[genome_tab$callular_mechanism == 'functional_apomixis',]
table(aggregate(ploidy ~ species, data=mitotic_parth, FUN=function(x){x[1]})$ploidy)

meiotic_parth <- genome_tab[!genome_tab$callular_mechanism %in% c('functional_apomixis', 'unknown'),]
table(aggregate(ploidy ~ species, data=meiotic_parth, FUN=function(x){x[1]})$ploidy)

hybrid_orig <- aggregate(hybrid_origin ~ species, data=genome_tab, FUN=function(x){x[1]})

# # genomes
cat("Number of reanalyzed species: ", number_of_species, "\n")

# genome sizes
cat("Range of haploid sizes [Mbp]:", paste(range(c(genome_tab[,'assembly_size[M]'], genome_tab[,'haploid_length[M]']), na.rm=T)), "\n")
cat("Haploid number of chromosomes: ", paste(range(literature_data$chromosomes_haploid, na.rm = T)), "\n")


# ploidy
cat("Number of diploid species: ", length(unique(genome_tab$species[diploids])), "\n")
cat("Number of polyploid species: ", length(unique(genome_tab$species[otherploids])), "\n")

cat("Max missing conserved genes in Meloidogyne: ", max(genome_tab[substr(genome_tab$species, 1, 11) == "Meloidogyne",'missing']), "\n")

cat("Max missing conserved genes in Meloidogyne: ", max(genome_tab[substr(genome_tab$species, 1, 11) == "Meloidogyne",'missing']), "\n")
cat("HGT in A. vaga: ", literature_data[literature_data$code == "Avag1", "HGT"], "\n")
cat("Heterozygosity of D. pachys: ", genome_tab[genome_tab$code == "Dpac1", 'heterozygosity'], "\n")

sum(table(genome_tab$reproduction_mode))
table(genome_tab[genome_tab$ploidy > 2,'reproduction_mode'])
sum(is.na(genome_tab$hybrid_origin))
genome_tab[genome_tab$ploidy > 2,c('species', 'hybrid_origin')]