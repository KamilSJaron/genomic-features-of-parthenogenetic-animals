# This function is taken form AsexStats package
#

filter_thetas <-
function (sp_data, min_cov = 0.5, filt_cov = T, window_size = 999,
    filt_window_size = T)
{
    if (filt_cov) {
        sp_data <- sp_data[sp_data$coverage > median(sp_data$coverage) *
            min_cov, ]
        sp_data <- sp_data[sp_data$coverage < (median(sp_data$coverage) *
            (1 + min_cov)), ]
    }
    if (filt_window_size) {
        sp_data <- sp_data[(sp_data$end - sp_data$start) >= window_size,
            ]
    }
    sp_data <- sp_data[!is.na(sp_data$theta_MLE), ]
    return(sp_data)
}
