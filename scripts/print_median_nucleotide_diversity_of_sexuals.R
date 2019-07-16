nt_div <- read.table('tables/nt_diversity.tsv', header = T, stringsAsFactors = F, sep = '\t')

get_population_mean <- function(entry){
    mean(as.numeric(unlist(strsplit(entry, ';'))))
}

nt_div <- nt_div[nt_div[,3] %in% c("Arthropoda", "Chordata", "Porifera", "Nematoda", "Mollusca", "Echinodermata"),]
nt_div <- as.vector(sapply(nt_div[,6], get_population_mean))
# print median of species (where value of a species is a mean of population estimates)
cat(round(quantile(nt_div, c(0.025, 0.975)), 2))
cat("\n")

# mean(sort(nt_div) < 0.03)
# [1] 0.03587444
# mean(sort(nt_div) < 0.53)
# [1] 0.4125561
