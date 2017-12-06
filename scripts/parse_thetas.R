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
plot_thetas <- function(i, title = '', col = 1, add = F){
    if (add){
        hist(log_thetas[[i]], breaks = 40, col = col, freq = F, add = T)
    } else {
        hist(log_thetas[[i]], breaks = 40, col = col, freq = F, xaxt="n",
             xlab = theta_label, main = title, ylim = c(0, 2.4), cex.lab = 1.3)
             at.x <- outer(1, 10^(c(-8, -3, -2, -1)))
             lab.x <- log10(at.x)
             axis(1, at=lab.x, labels=at.x, las=1, cex.lab = 1.3)
    }
}

png('figures/species_heterozygosity.png')
par(mfcol=c(3,4))
for(i in 1:length(theta_estimate_files)){
    plot_thetas(i, paste('s :', samples[i], ' r : ', references[i]))
}
dev.off()

col_1 = rgb(0, 0, 255, max = 255, alpha = 125)
col_2 = rgb(255, 0, 0, max = 255, alpha = 125)
# Mflo1
png('figures/thetas_Mflo.png')
    plot_thetas(8, add = F, col = col_1, title = 'Meloidogyne floridensis')
    plot_thetas(9, add = T, col = col_2)
    legend('topright', bty = 'n', col = c(col_1, col_2), legend = paste('s: ', samples[8:9], ' r: ', references[8:9]), pch = 20)
dev.off()

# M. incognita
png('figures/thetas_Minc.png')
    plot_thetas(10, add = F, col = col_1, title = 'Meloidogyne incognita')
    plot_thetas(11, add = T, col = col_2)
    legend('topright', bty = 'n', col = c(col_1, col_2), legend = paste('s: ', samples[10:11], ' r: ', references[10:11]), pch = 20)
dev.off()

# weighted.median(x, w, na.rm = TRUE)
library(matrixStats)

get_weighted_median <- function(x){
    wights <- x$end - x$start
    weightedMedian(x$theta_MLE, na.rm = T)
}
meadian_thetas <- unlist(lapply(atlas_data, get_weighted_median))
png('figures/median_theta.png', width=30, height=10)
    # MAKE median weighted by window sizes
    barplot(meadian_thetas, ylab = expression(paste("median " , theta)), names.arg = samples)
dev.off()
