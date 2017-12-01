source('scripts/filter_thetas.R')

theta_estimate_files <- list.files(path = ".", pattern = "theta_estimates.txt", recursive = T, include.dirs = T)

decomposed_file_names <- strsplit(theta_estimate_files, '/')
samples <- unlist(lapply(decomposed_file_names, function(x){x[2]}))
filenames <- unlist(lapply(decomposed_file_names, function(x){x[1]}))
references <- unlist(lapply(strsplit(filenames, split = '_'), function(x){ x[1] } ))

atlas_data <- list()
log_thetas <- list()

for(i in 1:length(theta_estimate_files)){
    file <- theta_estimate_files[i]
    atlas_data[[i]] <- read.table(file, header = T)
    log_thetas[[i]] <- log10(filter_thetas(atlas_data[[i]])[,'theta_MLE'])
}