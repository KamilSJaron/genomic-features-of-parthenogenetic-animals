source('scripts/filter_thetas.R')

theta_estimate_files <- list.files(path = ".", pattern = "theta_estimates.txt", recursive = T, include.dirs = T)

decomposed_file_names <- strsplit(theta_estimate_files, '/')
samples <- unlist(lapply(decomposed_file_names, function(x){x[2]}))
filenames <- unlist(lapply(decomposed_file_names, function(x){x[3]}))
references <- unlist(lapply(strsplit(filenames, split = '_'), function(x){ x[1] } ))

atlas_data <- list()
log_thetas <- list()

for(i in 1:length(theta_estimate_files)){
    file <- theta_estimate_files[i]
    atlas_data[[i]] <- read.table(file, header = T)
    log_thetas[[i]] <- log10(filter_thetas(atlas_data[[i]])[,'theta_MLE'])
}

theta_label <- expression(paste("[" , log[10], " " , theta, ']'))
plot_thetas <- function(i, title = ''){
    hist(log_thetas[[i]], breaks = 40, col = 1, freq = F, xaxt="n",
         xlab = theta_label, main = title, ylim = c(0, 2.4), cex.lab = 1.3)
    at.x <- outer(1, 10^(c(-8, -3, -2, -1)))
    lab.x <- log10(at.x)
    axis(1, at=lab.x, labels=at.x, las=1, cex.lab = 1.3)
}

png('figures/species_heterozygosity.png')
par(mfcol=c(3,4))
for(i in 1:length(theta_estimate_files)){
    plot_thetas(i, paste('s :', samples[i], ' r : ', references[i]))
}
dev.off()

