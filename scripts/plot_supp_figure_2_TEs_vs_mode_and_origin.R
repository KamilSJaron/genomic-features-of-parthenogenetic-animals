############
# Get data #
############

source('scripts/R_functions/load_genome_table.R')

genome_tab <- load_genome_table(c(1:18))
genome_tab <- genome_tab[nrow(genome_tab):1,]

hyb_origins <- c("no", "unknown", "yes")
repr_modes  <- c("gamete_duplication", "terminal_fusion", "central_fusion", "unknown_automixis",
                 "unknown", "functional_apomixis")

###
# misc
###
source('scripts/R_functions/load_cellular_mechanism_palette.R')
pal <- load_cellular_mechanism_palette()

plot_corpus <- function(.var, .ylab, .legend = F){
    plot(NULL, bty = 'n', axes=FALSE, xlim = c(0.65, 3.35), ylim = c(0,max(genome_tab[,.var], na.rm = T)),
        xlab = "Hybrid origin", ylab = .ylab,
        pch = NA, col = rgb(0.878,0.878,0.878), cex.lab = 1.3
    )

    axis(2, las=1)
    axis(1, labels = hyb_origins, at = 1:3, tick = F, line = F, cex.axis = 1.4)
    rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col=rgb(0.92,0.92,0.92))

    boxplot_col <- rgb(0.66, 0.66, 0.66, alpha=0.5)
    known_hybrid_orig_tab <- genome_tab[genome_tab$hybrid_origin != 'unknown',]
    boxplot(known_hybrid_orig_tab[,.var] ~ known_hybrid_orig_tab$hybrid_origin, bty = 'n', axes=FALSE,
          pch = NA, col = boxplot_col, cex.lab = 1.3, add = T
    )

    # grid(lwd = 1.2, col = 'gray', nx=NA, ny=NULL)
    grid(lwd = 1, col = 1, nx=3, ny=NULL)
    abline(lwd = 0.1, h = seq(0, max(genome_tab[,.var], na.rm = T), by = 5)[(1:28) %% 4 != 1])

    if(.legend){
        legend('topleft', bty = 'n',
               c("gamete duplication", "terminal fusion", "central fusion", "unknown meiosis", "unknown", "functional mitosis"),
               col = pal, pch = c(rep(20,6)), cex = 1)
    }
}

plot_ploits <- function(hyb_origin = "no", .var = "repeats"){
    # repr_modes
    subset <- genome_tab[genome_tab$hybrid_origin == hyb_origin, c(.var, "ploidy", "callular_mechanism")]
    subset <- subset[order(subset$callular_mechanism),]
    subset <- subset[!is.na(subset[,.var]),]

    # if(hyb_origin == "unknown"){
    #     subset <- subset[c(1,4,2,3,5:8),]
    # }

    at <- which(hyb_origin == hyb_origins)
    misplacement <- seq(from = -0.35, to = 0.35,
                        length = nrow(subset) + 2)[2:(nrow(subset)+1)]

    points(at + misplacement, subset[,.var],
           col = pal[subset$callular_mechanism],
           pch = 19, cex = 1.3)
    points(at + misplacement, subset[,.var],
           pch = 21, cex = 1.3)
    text(at + misplacement, subset[,.var], substr(rownames(subset), 1, 4), pos = 1)
}

#####
# REPEATS
####


tiff("figures/Supp_fig2_rep_TE_patterns.tiff",
     width = 12, height = 6, units = 'in', res = 90)
# pdf('figures/Supp_fig2_rep_TE_patterns.pdf', width = 16, height = 8)

# tiff("figures/Supp_fig2a_repetitions.tiff",
#      width = 8, height = 8, units = 'in', res = 90)
# # png('figures/Supp_fig2a_repetitions.png')

# plot_tab <- genome_tab[!is.na(genome_tab$repeats),]
# # par(mfrow = c(1,2))
# plot_corpus("repeats", "Repetitions [%]", T)
# title("Overall repetitive content")
# plot_ploits("no", "repeats")
# plot_ploits("unknown", "repeats")
# plot_ploits("yes", "repeats")

# dev.off()

#####
# TEs
####
plot_tab <- genome_tab[!is.na(genome_tab$TEs),]

tiff("figures/Supp_fig2b_TEs.tiff",
     width = 8, height = 8, units = 'in', res = 90)
# # png('figures/Supp_fig2b_TEs.png')

plot_corpus("TEs", "TEs [%]", T)
title("Transposable elements")
plot_ploits("no", "TEs")
plot_ploits("unknown", "TEs")
plot_ploits("yes", "TEs")

dev.off()
