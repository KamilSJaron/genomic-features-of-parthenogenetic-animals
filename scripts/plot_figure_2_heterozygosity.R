#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

require('gsheet')
require('RColorBrewer')
library('extrafont')
library('shape')

loadfonts(quiet=TRUE)

# if flag --presentation is specified; a figre for presentation is generated
presentation <- F
if ( length(args) == 1 ) {
  if ( args[1] == "--presentation" ){
    presentation = T
  }
}

############
# Get data #
############

# literature_data <- read.csv(text = gsheet2text("https://docs.google.com/spreadsheets/d/1T4BHQxzMGMlWiNJ7G9OLPXzdNBCqMw0mSO7C_ggGEbc/edit?usp=sharing", format='csv'),
#              stringsAsFactors = F, skip = 1, header = T, check.names = F)
# reported_heterozygosity <- literature_data[,c("code","haplotype_divergence[%]")]

genome_tab <- read.table('tables/genome_table.tsv',
                         header = T, stringsAsFactors = F, skip = 1, check.names = F)
rownames(genome_tab) <- genome_tab$code

# genome_tab['Rmac1','heterozygosity'] <- 0.123 # these estimates should be added to the big genome table, not here
# genome_tab['Rmag1','heterozygosity'] <- 0.492
# genome_tab['Dcor1','callular_mechanism'] <- "unknown_automixis"
# genome_tab[4:5,15:16] <- NA # delete rotifers

genome_tab$callular_mechanism[is.na(genome_tab$callular_mechanism)] <- "unknown"
genome_tab$hybrid_origin[is.na(genome_tab$hybrid_origin)] <- "unknown"

hyb_origins <- c("no", "unknown", "yes")
repr_modes  <- c("gamete_duplication", "terminal_fusion", "central_fusion", "unknown_automixis",
                 "unknown", "functional_apomixis")

# ordering
genome_tab$hybrid_origin <- ordered(genome_tab$hybrid_origin, levels=hyb_origins )
genome_tab$callular_mechanism <- ordered(genome_tab$callular_mechanism, levels=repr_modes )
genome_tab['Avag1','heterozygosity'] <- NA
genome_tab <- genome_tab[!is.na(genome_tab$heterozygosity),]

########
# plot #
########

source('scripts/R_functions/load_cellular_mechanism_palette.R')
pal <- load_cellular_mechanism_palette()

plot_corpus <- function(){
    plot(NULL, bty = 'n', axes=FALSE, xlim = c(0.65, 3.35), ylim = c(0,12.4),
        xlab = "Hybrid origin", ylab = "Heterozygosity [%]",
        pch = NA, col = rgb(0.878,0.878,0.878), cex.lab = 1.4
    )

    axis(2, las=1)
    axis(1, labels = hyb_origins, at = 1:3, tick = F, line = F, cex.axis = 1.4)
    rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col=rgb(0.92,0.92,0.92))

    ellipse_col <- rgb(0.66, 0.66, 0.66, alpha=0.5)
    filledellipse(rx1 = 3.5, ry1 = 0.21, col = ellipse_col, angle = 90.8, dr = 0.1, mid = c(1.97, 9.3))
    filledellipse(rx1 = 3.8, ry1 = 0.2, col = ellipse_col, angle = 87.5, dr = 0.1, mid = c(3.13, 6))
    filledellipse(rx1 = 0.2, ry1 = 0.9, col = ellipse_col, angle = 0, dr = 0.1, mid = c(2.67, 4.86))

    boxplot_col <- rgb(0.66, 0.66, 0.66, alpha=0.5)
    known_hybrid_orig_tab <- genome_tab[genome_tab$hybrid_origin != 'unknown',]
    with(known_hybrid_orig_tab,
             boxplot(heterozygosity ~ hybrid_origin, bty = 'n', axes=FALSE,
                      xlab = "Hybrid origin", ylab = "Heterozygosity [%]",
                      pch = NA, col = boxplot_col, cex.lab = 1.4, add = T
             )
    )

    # grid(lwd = 1.2, col = 'gray', nx=NA, ny=NULL)
    grid(lwd = 2, col = 1, nx=3, ny=NA)
    grid(lwd = 0.5, col = 1, nx=NA, ny=NULL, lty = 3)
    # abline(lwd = 0.1, h = seq(0, 14, by = 0.5)[(1:28) %% 4 != 1])

    legend('topleft', bty = 'n',
           c("gamete duplication", "terminal fusion", "central fusion", "unknown meiosis", "unknown", "functional mitosis"),
           col = pal, pch = c(rep(20,6)), cex = 1.2)
    # legend('topleft', bty = 'n',
    #        c("gamete duplication", "terminal fusion", "central fusion", "unknown automixis", "unknown", "functional apomixis",
    #          '','diploid','triploid', 'tetraploid'),
    #        col = c(pal, NA, 1, 1, 1), pch = c(rep(20,6), NA, 19, 17, 15), cex = 1.2)
}

plot_ploits <- function(hyb_origin = "no"){
    # repr_modes
    subset <- genome_tab[genome_tab$hybrid_origin == hyb_origin, c("heterozygosity", "ploidy", "callular_mechanism")]
    subset <- subset[order(subset$callular_mechanism),]
    subset <- subset[!is.na(subset$heterozygosity),]

    if(hyb_origin == "unknown"){
        subset <- subset[c(1,4,2,3,5:8),]
    }

    at <- which(hyb_origin == hyb_origins)
    misplacement <- seq(from = -0.35, to = 0.35, length = nrow(subset) + 2)[2:(nrow(subset)+1)]
    # symbols <- c(NA, 19, 17, 15)
    # conture_symbols <- c(NA, 21, 24, 22)

    points(at + misplacement, subset$heterozygosity,
           col = pal[subset$callular_mechanism],
           pch = 19, cex = 1.3) #symbols[subset$ploidy]
    points(at + misplacement, subset$heterozygosity,
           pch = 21, cex = 1.3) #conture_symbols[subset$ploidy]
    # for(line_number in 1:nrow(subset)){
    #     lines(rep(at,2) + misplacement[line_number], subset[line_number, 1:2])
    # }
}

pdf('figures/fig2a_heterozygosity.pdf', width=8, height=8)
# png('figures/fig2_heterozygosity.png')

# 'c(bottom, left, top, right)'
par(mar = c(4.5, 4.5, 2, 1) + 0.1)
plot_corpus()

for ( ho in hyb_origins ){
    plot_ploits(ho)
}

dev.off()














