tab_file <- 'tables/genome_table.tsv'
genome_tab <- read.table(tab_file, header = T, stringsAsFactors = F, skip = 1, check.names = F)

# get unique species
genome_tab <- genome_tab[!is.na(genome_tab$complete),]

# arthropods and vertebrates
range(rowSums(genome_tab[genome_tab$code %in% c('Pfor1', 'Lcla1', 'Tpre1', 'Obir1', 'Aruf1', 'Fcan1', 'Pvir1'), c('complete', 'fragmented')]))

# rotifers
range(genome_tab[genome_tab$code %in% c('Avag1', 'Aric1', 'Rmac1', 'Rmag1'), c('missing')])

# nematodes
nematodes <- genome_tab[genome_tab$complete < 80 | genome_tab$code %in% c('Dcor1','Dpac1'),c('complete', 'fragmented', 'duplicated', 'missing')]