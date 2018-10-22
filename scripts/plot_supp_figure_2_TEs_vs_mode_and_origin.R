############
# Get data #
############

genome_tab <- read.table('tables/genome_table.tsv',
                         header = T, stringsAsFactors = F, skip = 1, check.names = F)
rownames(genome_tab) <- genome_tab$code
# take only relevant columns and reverse the order of rows
genome_tab <- genome_tab[nrow(genome_tab):1,c(1,3:5,15:18)]
# kick out those that have no
genome_tab <- genome_tab[!(is.na(genome_tab$TEs) & is.na(genome_tab$repeats)),]

genome_tab$reproduction_mode[is.na(genome_tab$reproduction_mode)] <- "unknown"
genome_tab$hybrid_origin[is.na(genome_tab$hybrid_origin)] <- "unknown"
hyb_origins <- c("no", "unknown", "yes")
repr_modes  <- c("gamete_duplication", "terminal_fusion", "central_fusion", "unknown_automixis",
                 "unknown", "functional_apomixis")

# ordering
genome_tab$hybrid_origin <- ordered(genome_tab$hybrid_origin, levels=hyb_origins )
genome_tab$reproduction_mode <- ordered(genome_tab$reproduction_mode, levels=repr_modes )

###
# misc
###

pal <- c(rgb(0.9411, 0.8941, 0.2588), # yellow - gamete duplication
         rgb(0.9019, 0.6235, 0     ), # orange - terminal fusion
         rgb(0.8352, 0.3686, 0     ), # vermillion - central fusion
         rgb(0.8,    0.4745, 0.6549), # pink - automixis unknown
         rgb(0,      3/5,    2/5     ), # blue - apomixis
         rgb(0,      0.4470, 0.6980)) # green - unknown

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
               c("gamete duplication", "terminal fusion", "central fusion", "unknown automixis", "unknown", "functional apomixis"),
               col = pal, pch = c(rep(20,6)), cex = 1)
    }
}

plot_ploits <- function(hyb_origin = "no", .var = "repeats"){
    # repr_modes
    subset <- genome_tab[genome_tab$hybrid_origin == hyb_origin, c(.var, "ploidy", "reproduction_mode")]
    subset <- subset[order(subset$reproduction_mode),]
    subset <- subset[!is.na(subset[,.var]),]

    # if(hyb_origin == "unknown"){
    #     subset <- subset[c(1,4,2,3,5:8),]
    # }

    at <- which(hyb_origin == hyb_origins)
    misplacement <- seq(from = -0.35, to = 0.35,
                        length = nrow(subset) + 2)[2:(nrow(subset)+1)]

    points(at + misplacement, subset[,.var],
           col = pal[subset$reproduction_mode],
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

plot_tab <- genome_tab[!is.na(genome_tab$repeats),]
par(mfrow = c(1,2))
plot_corpus("repeats", "Repetitions [%]", T)
title("Overall repetitive content")
plot_ploits("no", "repeats")
plot_ploits("unknown", "repeats")
plot_ploits("yes", "repeats")

# dev.off()

#####
# TEs
####
plot_tab <- genome_tab[!is.na(genome_tab$TEs),]

# tiff("figures/Supp_fig2b_TEs.tiff",
#      width = 8, height = 8, units = 'in', res = 90)
# # png('figures/Supp_fig2b_TEs.png')

plot_corpus("TEs", "TEs [%]")
title("Transposable elements")
plot_ploits("no", "TEs")
plot_ploits("unknown", "TEs")
plot_ploits("yes", "TEs")

dev.off()
