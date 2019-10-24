#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

require('gsheet')
require('RColorBrewer')
library('extrafont')
library('shape')
library('plotrix')

source('scripts/R_functions/output_parsers.R')

loadfonts(quiet=TRUE)

# if flag --presentation is specified; a figre for presentation is generated
presentation <- T
split_axis <- T
homoeolog <- T
rm_boxes <- T

############
# Get data #
############

# literature_data <- read.csv(text = gsheet2text("https://docs.google.com/spreadsheets/d/1T4BHQxzMGMlWiNJ7G9OLPXzdNBCqMw0mSO7C_ggGEbc/edit?usp=sharing", format='csv'),
#              stringsAsFactors = F, skip = 1, header = T, check.names = F)
# reported_heterozygosity <- literature_data[,c("code","haplotype_divergence[%]")]

genome_tab <- read.table('tables/genome_table.tsv',
                         header = T, stringsAsFactors = F, skip = 1, check.names = F, sep = '\t')
rownames(genome_tab) <- genome_tab$code

genome_tab[c('Tdi1', 'Tsi1', 'Tms1', 'Tte1', 'Tge1'),'heterozygosity'] <- c(0.107, 0.0823, 0.0859, 0.0878, 0.103)
genome_tab[c('Tdi1', 'Tsi1', 'Tms1', 'Tte1', 'Tge1'),'callular_mechanism'] <- "functional_apomixis"
    # Tdi 0.107
    # Tsi 0.0823
    # Tms 0.0859
    # Tte 0.0878
    # Tge 0.103

# genome_tab['Rmac1','heterozygosity'] <- 0.123 # these estimates should be added to the big genome table, not here
# genome_tab['Rmag1','heterozygosity'] <- 0.492
# genome_tab['Dcor1','callular_mechanism'] <- "unknown_automixis"
# genome_tab[4:5,15:16] <- NA # delete rotifers

genome_tab$callular_mechanism[is.na(genome_tab$callular_mechanism)] <- "unknown"

repr_modes  <- c("gamete_duplication", "terminal_fusion", "central_fusion", "unknown_automixis",
                 "unknown", "functional_apomixis")

rotifers <- c('Rmac1', 'Rmag1', 'Avag1', 'Aric1')
# ./data/*/genomescope_diploid
# medians from email 12.7.2019
# Summary values for % divergence between ohnologous pairs of genes are:
#   A. ricciae: mean = 33.5254 ± 8.64032 SD (median = 33.2108), N = 10533
#   A. vaga: mean = 28.5242 ± 8.51263 SD (median = 27.8912), N = 8577
#   R. macrura: mean = 29.2332 ± 8.0746 SD (median = 28.9116), N = 1255
#   R. magnacalcarata: mean = 27.5818 ± 8.15129 SD (median = 27.0777), N = 1400
rotifers_ohno <- c(28.91, 27.08, 27.89, 33.21)

diploid_rot_genomescope_files <- paste0('./data/', rotifers, '/genomescope_diploid/', rotifers, '_summary.txt')
rotifers_homo <- sapply(diploid_rot_genomescope_files, function(x) { parse_genomescope_summary(x)[3] })
rotifers_total <- rotifers_ohno + rotifers_homo
genome_tab[rotifers,'heterozygosity'] <- 50
genome_tab[rotifers,'callular_mechanism'] <- 'functional_apomixis'

g_from <- 10.3
g_to <- 26
values_to_shift <- which(genome_tab$heterozygosity > g_to)
# genome_tab[values_to_shift,'heterozygosity'] <- genome_tab[values_to_shift,'heterozygosity']

genome_tab <- genome_tab[!is.na(genome_tab$heterozygosity),]

# ordering
genome_tab$callular_mechanism <- ordered(genome_tab$callular_mechanism, levels=repr_modes )


########
# plot #
########

source('scripts/R_functions/load_cellular_mechanism_palette.R')
pal <- load_cellular_mechanism_palette()
ellipse_col <- rgb(0.66, 0.66, 0.66, alpha=0.5)

plot_corpus <- function(presentation = F){
    general_cex <- ifelse(presentation, 1.6, 1.4)
    ymax <- ifelse(homoeolog, 33.2, 39)
    gap.plot(100, c(-5), gap=c(g_from, g_to), xlim = c(0.65, 3.35), ylim = c(0, ymax),
        xlab = " ", ylab = "Heterozygosity [%]", cex.lab = general_cex)
    rect(par("usr")[1], par("usr")[3], par("usr")[2], g_from, col=rgb(0.92,0.92,0.92))
    rect(par("usr")[1], g_from*(1+0.02), par("usr")[2], par("usr")[4], col=rgb(0.92,0.92,0.92))
    axis.break(2, g_from, breakcol="snow", style="gap")
    axis.break(2, g_from*(1+0.02), breakcol="black", style="slash")
    axis.break(4, g_from*(1+0.02), breakcol="black", style="slash")
    axis(2, seq(0, g_from, by = 2), col = NA, col.ticks = 1)
    axis(2, seq(12, ymax - 15, by = 2), col = NA, col.ticks = 1, labels = seq(27, ymax, by = 2))

    # axis(1, labels = hyb_origins, at = 1:3, tick = F, line = F, cex.axis = general_cex)

    filledellipse(rx1 = 13, ry1 = 0.18, col = ellipse_col, angle = 89.7, dr = 0.1, mid = c(3.175, 9.4))
    # Meloidogyne
    filledellipse(rx1 = 2.45, ry1 = 0.2, col = ellipse_col, angle = 88.5, dr = 0.1, mid = c(2.78, 7.23))
    # diploscapter
    filledellipse(rx1 = 0.16, ry1 = 1.5, col = ellipse_col, angle = 0, dr = 0.1, mid = c(1.52, 4.6))

    legend('topleft', bty = 'n',
           c("gamete duplication",
             "terminal fusion",
             "central fusion",
             "unknown meiosis",
             "unknown",
             "functional mitosis"),
           col = pal, pch = c(rep(20,6)), cex = ifelse(presentation, 1.5, 1.2))
}

point_size <- ifelse(presentation, 1.5, 1.3)

fig_file <- 'figures/presentation/fig2_heterozygosity_by_cellular_mode.pdf'
print(paste("Saving...", fig_file))
fig_width <- ifelse(presentation, 12, 10)

pdf(fig_file, width = fig_width, height = 8)
# png('figures/fig2_heterozygosity.png')

# 'c(bottom, left, top, right)'
par(mar = c(4.5, 4.5, 2, 1) + 0.1)
plot_corpus(presentation)

# repr_modes
genome_tab <- genome_tab[, c("heterozygosity", "ploidy", "callular_mechanism")]
genome_tab <- genome_tab[!is.na(genome_tab$heterozygosity),]
genome_tab <- genome_tab[order(genome_tab$heterozygosity),]
genome_tab <- genome_tab[order(genome_tab$callular_mechanism),]

placement <- seq(from = 0.7, to = 3.3, length = nrow(genome_tab))
# # symbols <- c(NA, 19, 17, 15)
# # conture_symbols <- c(NA, 21, 24, 22)
#
composite <- which(genome_tab$ploidy > 2)
# # if ( hyb_origin == "yes" ){
points(placement, genome_tab$heterozygosity,
       bg = pal[genome_tab$callular_mechanism],
       pch = 24, cex = point_size)
points((placement)[composite], genome_tab$heterozygosity[composite],
       bg = pal[genome_tab$callular_mechanism][composite],
       pch = 25, cex = point_size)

xpos <- placement[row.names(genome_tab) %in% c('Avag1', 'Aric1', 'Rmac1', 'Rmag1')]
# symbols <- c(NA, 19, 17, 15)
# conture_symbols <- c(NA, 21, 24, 22)
rotifers_to_plot <- rotifers_ohno - (g_to - g_from)

points(xpos, rotifers_to_plot,
       bg = pal[6],
       pch = 24, cex = point_size) #symbols[subset$ploidy]
# points(xpos, rotifers_to_plot,
#        pch = 24, cex = point_size) #conture_symbols[subset$ploidy]

points(xpos, rotifers_homo,
       bg = pal[6],
       pch = 25, cex = point_size) #symbols[subset$ploidy]
# points(xpos, rotifers_homo,
#        pch = 25, cex = point_size) #conture_symbols[subset$ploidy]

for ( i in 1:4){
    lines(c(xpos[i], xpos[i]), c(rotifers_homo[i], rotifers_to_plot[i]), lty = 5, lwd = 0.6)
}

legend(
    2.13, 18.5, bty = 'n', cex = 1.5,
    paste(c('homoeolog', 'composite', 'allelic'), 'heterozygosity'),
    pch = c(24, 24, 25)
)
legend(
    2.13, 18.5, bty = 'n', cex = 1.5,
    rep(' ', 3),
    pch = c(NA, 25, NA)
)

dev.off()


