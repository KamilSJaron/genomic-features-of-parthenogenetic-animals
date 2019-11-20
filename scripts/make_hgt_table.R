# pull data from Reuben's google doc and turn it into supplementary a table

library('gsheet')

hgt_data <- read.csv(text = gsheet2text("https://docs.google.com/spreadsheets/d/1vxEQ51UdlunRDR9Nr9eeDGI8xSHGd-Hqa9fq8C3fzv0/edit?usp=sharing", format='csv'),
                     stringsAsFactors = F, header = T, check.names = F)

columns_to_keep <- c('id', 'num_proteins', 'uniref_hits', 'uniref_hits_perc',
                     'kphylum', 'kphylum%', 'kphylum+linked', 'kphylum+linked%')

# keep only the complete data
hgt_data <- hgt_data[hgt_data$status == 'complete',]

sexuals <- grepl("^s", hgt_data$id)
hgt_data_sexuals <- hgt_data[sexuals, ]
hgt_data_asexuals <- hgt_data[!sexuals, ]

length(unique(hgt_data_asexuals$species)) - 2
# - 2 is for three Panagrolaimus strains that are considered as one species in this study
# 23 species