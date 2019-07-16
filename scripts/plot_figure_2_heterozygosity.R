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
presentation <- ifelse( "--presentation" %in% args, T, F)
excl_rotifers <- ifelse( "--excl_rotifers" %in% args, T, F)
roti_arrow <- ifelse( "--arrow" %in% args, T, F)
split_axis <- ifelse( "--split_axis" %in% args, T, F)
homoeolog <- ifelse( "--homoeolog" %in% args, T, F)

############
# Get data #
############

# literature_data <- read.csv(text = gsheet2text("https://docs.google.com/spreadsheets/d/1T4BHQxzMGMlWiNJ7G9OLPXzdNBCqMw0mSO7C_ggGEbc/edit?usp=sharing", format='csv'),
#              stringsAsFactors = F, skip = 1, header = T, check.names = F)
# reported_heterozygosity <- literature_data[,c("code","haplotype_divergence[%]")]

genome_tab <- read.table('tables/genome_table.tsv',
                         header = T, stringsAsFactors = F, skip = 1, check.names = F, sep = '\t')
rownames(genome_tab) <- genome_tab$code


if ( presentation ){
    presentation = T
    genome_tab[c('Tdi1', 'Tsi1', 'Tms1', 'Tte1', 'Tge1'),'heterozygosity'] <- c(0.107, 0.0823, 0.0859, 0.0878, 0.103)
    genome_tab[c('Tdi1', 'Tsi1', 'Tms1', 'Tte1', 'Tge1'),'hybrid_origin'] <- "no"
    genome_tab[c('Tdi1', 'Tsi1', 'Tms1', 'Tte1', 'Tge1'),'callular_mechanism'] <- "functional_apomixis"
    # Tdi 0.107
    # Tsi 0.0823
    # Tms 0.0859
    # Tte 0.0878
    # Tge 0.103
}

# genome_tab['Rmac1','heterozygosity'] <- 0.123 # these estimates should be added to the big genome table, not here
# genome_tab['Rmag1','heterozygosity'] <- 0.492
# genome_tab['Dcor1','callular_mechanism'] <- "unknown_automixis"
# genome_tab[4:5,15:16] <- NA # delete rotifers

genome_tab$callular_mechanism[is.na(genome_tab$callular_mechanism)] <- "unknown"
genome_tab$hybrid_origin[is.na(genome_tab$hybrid_origin)] <- "unknown"

hyb_origins <- c("no", "unknown", "yes")
repr_modes  <- c("gamete_duplication", "terminal_fusion", "central_fusion", "unknown_automixis",
                 "unknown", "functional_apomixis")

rotifers <- c('Rmac1', 'Rmag1', 'Avag1', 'Aric1')
if ( excl_rotifers ) {
    genome_tab[rotifers, 'heterozygosity'] <- NA
} else {
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

    # whatever will be the homoeolog divergence
    if ( homoeolog ) {
        genome_tab[rotifers,'heterozygosity'] <- 50
    } else {
        genome_tab[rotifers,'heterozygosity'] <- rotifers_total
    }
}

if ( split_axis ){
    g_from <- 10.3
    g_to <- 26
    values_to_shift <- which(genome_tab$heterozygosity > g_to)
    genome_tab[values_to_shift,'heterozygosity'] <- genome_tab[values_to_shift,'heterozygosity']
    # print(genome_tab[rotifers,'heterozygosity'])
}

genome_tab <- genome_tab[!is.na(genome_tab$heterozygosity),]

# ordering
genome_tab$hybrid_origin <- ordered(genome_tab$hybrid_origin, levels=hyb_origins )
genome_tab$callular_mechanism <- ordered(genome_tab$callular_mechanism, levels=repr_modes )

########
# plot #
########

source('scripts/R_functions/load_cellular_mechanism_palette.R')
pal <- load_cellular_mechanism_palette()
ellipse_col <- rgb(0.66, 0.66, 0.66, alpha=0.5)

plot_corpus <- function(presentation = F){
    general_cex <- ifelse(presentation, 1.6, 1.4)
    if ( split_axis ){
        ymax <- ifelse(homoeolog, 33.2, 39)
        gap.plot(100, c(-5), gap=c(g_from, g_to), xlim = c(0.65, 3.35), ylim = c(0, ymax),
            xlab = "Hybrid origin", ylab = "Divergence [%]", cex.lab = general_cex)
        rect(par("usr")[1], par("usr")[3], par("usr")[2], g_from, col=rgb(0.92,0.92,0.92))
        rect(par("usr")[1], g_from*(1+0.02), par("usr")[2], par("usr")[4], col=rgb(0.92,0.92,0.92))
        axis.break(2, g_from, breakcol="snow", style="gap")
        axis.break(2, g_from*(1+0.02), breakcol="black", style="slash")
        axis.break(4, g_from*(1+0.02), breakcol="black", style="slash")
        axis(2, seq(0, g_from, by = 2), col = NA, col.ticks = 1)
        axis(2, seq(12, ymax - 15, by = 2), col = NA, col.ticks = 1, labels = seq(27, ymax, by = 2))
    } else {
        ymax <- ifelse(roti_arrow, 9.75, 12.4)
        plot(NULL, bty = 'n', axes=FALSE, xlim = c(0.65, 3.35), ylim = c(0, ymax),
            xlab = "Hybrid origin", ylab = "Divergence [%]",
            pch = NA, col = rgb(0.878,0.878,0.878), cex.lab = general_cex)
            axis(2, las = 1, lwd = 0)
        rect(par("usr")[1], par("usr")[3], par("usr")[2], par("usr")[4], col=rgb(0.92,0.92,0.92))
        axis(2, las=1)
    }

    axis(1, labels = hyb_origins, at = 1:3, tick = F, line = F, cex.axis = general_cex)

    if ( presentation) {
      filledellipse(rx1 = 3.5, ry1 = 0.21, col = ellipse_col, angle = 90.8, dr = 0.1, mid = c(1.97, 9.3))
      filledellipse(rx1 = 3.8, ry1 = 0.2, col = ellipse_col, angle = 87.5, dr = 0.1, mid = c(3.13, 6))
      filledellipse(rx1 = 0.15, ry1 = 0.93, col = ellipse_col, angle = 0, dr = 0.1, mid = c(2.68, 4.87))
    } else {
      if ( !excl_rotifers & ! split_axis & ! roti_arrow ){
        # rotifers
        filledellipse(rx1 = 3.5, ry1 = 0.21, col = ellipse_col, angle = 90.8, dr = 0.1, mid = c(1.97, 9.3))
      }

      if ( ! homoeolog){
          # Meloidogyne
          filledellipse(rx1 = 3.8, ry1 = 0.2, col = ellipse_col, angle = 87.5, dr = 0.1, mid = c(3.13, 6))
          # diploscapter
          filledellipse(rx1 = 0.175, ry1 = 0.9, col = ellipse_col, angle = 0, dr = 0.1, mid = c(2.67, 4.86))
      }
    }
    if ( split_axis ){
        if ( homoeolog ){
            filledellipse(rx1 = 3.6, ry1 = 0.28, col = ellipse_col, angle = 90, dr = 0.1, mid = c(2.05, 14.55))
            # Meloidogyne
            filledellipse(rx1 = 3.8, ry1 = 0.24, col = ellipse_col, angle = 86.5, dr = 0.1, mid = c(3.17, 6))
            # diploscapter
            filledellipse(rx1 = 0.16, ry1 = 1.3, col = ellipse_col, angle = 0, dr = 0.1, mid = c(2.67, 4.6))
        } else {
            filledellipse(rx1 = 6.5, ry1 = 0.25, col = ellipse_col, angle = 90, dr = 0.1, mid = c(2.06, 17.5))
        }
    }

    boxplot_col <- rgb(0.66, 0.66, 0.66, alpha=0.5)
    known_hybrid_orig_tab <- genome_tab[genome_tab$hybrid_origin != 'unknown',]
    with(known_hybrid_orig_tab,
             boxplot(heterozygosity ~ hybrid_origin, bty = 'n', axes=FALSE,
                      xlab = "Hybrid origin", ylab = "Divergence [%]",
                      pch = NA, col = boxplot_col, cex.lab = general_cex, add = T
             )
    )

    # grid(lwd = 1.2, col = 'gray', nx=NA, ny=NULL)
    grid(lwd = 2, col = 1, nx=3, ny=NA)
    if ( !split_axis ){
        grid(lwd = 0.5, col = 1, nx=NA, ny=NULL, lty = 3)
    }
    # abline(lwd = 0.1, h = seq(0, 14, by = 0.5)[(1:28) %% 4 != 1])

    legend('topleft', bty = 'n',
           c("gamete duplication",
             "terminal fusion",
             "central fusion",
             "unknown meiosis",
             "unknown",
             "functional mitosis"),
           col = pal, pch = c(rep(20,6)), cex = ifelse(presentation, 1.6, 1.2))
    # legend('topleft', bty = 'n',
    #        c("gamete duplication", "terminal fusion", "central fusion", "unknown automixis", "unknown", "functional apomixis",
    #          '','diploid','triploid', 'tetraploid'),
    #        col = c(pal, NA, 1, 1, 1), pch = c(rep(20,6), NA, 19, 17, 15), cex = 1.2)
}

point_size <- ifelse(presentation, 1.5, 1.3)

plot_ploits <- function(hyb_origin = "no", presentation = F){
    # repr_modes
    subset <- genome_tab[genome_tab$hybrid_origin == hyb_origin, c("heterozygosity", "ploidy", "callular_mechanism")]
    subset <- subset[order(subset$callular_mechanism),]
    subset <- subset[!is.na(subset$heterozygosity),]

    if(hyb_origin == "unknown" & !excl_rotifers){
        subset <- subset[c(1,2,5,3,4,6:9),]
    }

    at <- which(hyb_origin == hyb_origins)
    misplacement <- seq(from = -0.45, to = 0.45, length = nrow(subset) + 2)[2:(nrow(subset)+1)]
    # symbols <- c(NA, 19, 17, 15)
    # conture_symbols <- c(NA, 21, 24, 22)

    if ( homoeolog ){
        composite <- which(subset$ploidy > 2)
        points((at + misplacement)[composite], subset$heterozygosity[composite],
               bg = pal[subset$callular_mechanism][composite],
               pch = 24, cex = point_size)
        points(at + misplacement, subset$heterozygosity,
               bg = pal[subset$callular_mechanism],
               pch = 25, cex = point_size)
    } else {
        points(at + misplacement, subset$heterozygosity,
               col = pal[subset$callular_mechanism],
               pch = 19, cex = point_size) #symbols[subset$ploidy]
        points(at + misplacement, subset$heterozygosity,
               pch = 21, cex = point_size) #conture_symbols[subset$ploidy]
    }
    # for(line_number in 1:nrow(subset)){
    #     lines(rep(at,2) + misplacement[line_number], subset[line_number, 1:2])
    # }
}

fig_file <- paste0('figures/fig2a_heterozygosity',
                   ifelse(presentation, '_presentation', ''),
                   ifelse(excl_rotifers, '_excl_rotifers', ''),
                   ifelse(roti_arrow, '_roti_arrow', ''),
                   ifelse(split_axis, '_split_axis', ''),
                   ifelse(homoeolog, '_homoeolog', ''),
                   '.pdf')
print(paste("Saving...", fig_file))
fig_width <- ifelse(presentation, 12, 10)

pdf(fig_file, width = fig_width, height = 8)
# png('figures/fig2_heterozygosity.png')

# 'c(bottom, left, top, right)'
par(mar = c(4.5, 4.5, 2, 1) + 0.1)
plot_corpus(presentation)

for ( ho in hyb_origins ){
    plot_ploits(ho, presentation)
}


if ( homoeolog ){
    unknown_repr <- sum(genome_tab$hybrid_origin == "unknown")

    misplacement <- seq(from = -0.45, to = 0.45, length = unknown_repr + 2)[2:(unknown_repr+1)]
    print(misplacement)
    xpos <- 2 + misplacement[c(4,5,6,7)]
    # symbols <- c(NA, 19, 17, 15)
    # conture_symbols <- c(NA, 21, 24, 22)
    rotifers_to_plot <- rotifers_ohno - (g_to - g_from)

    points(xpos, rotifers_to_plot,
           bg = pal[c(5,5,6,6)],
           pch = 24, cex = point_size) #symbols[subset$ploidy]
    # points(xpos, rotifers_to_plot,
    #        pch = 24, cex = point_size) #conture_symbols[subset$ploidy]

    points(xpos, rotifers_homo,
           bg = pal[c(5,5,6,6)],
           pch = 25, cex = point_size) #symbols[subset$ploidy]
    # points(xpos, rotifers_homo,
    #        pch = 25, cex = point_size) #conture_symbols[subset$ploidy]

    for ( i in 1:4){
        lines(c(xpos[i], xpos[i]), c(rotifers_homo[i], rotifers_to_plot[i]), lty = 5, lwd = 0.6)
    }
}

if ( roti_arrow ){
    filledellipse(rx1 = 0.35, ry1 = 1.1, col = ellipse_col, angle = 0, dr = 0.1, mid = c(2, 9))
    arrows(1.8, 8.2, 1.8, 9.6,
        code = 2, col = par("fg"), lty = par("lty"), lwd = 3)
    text(rep(2.08, 4), seq(8.3, 9.5, len = 4), paste(rotifers, rotifers_total, "%"), cex = 1.2, las = 1)
}

if ( homoeolog ){
    legend(
        0.54, 13.5, bty = 'n', cex = ifelse(presentation, 1.6, 1.2),
        paste(c('homoeolog', 'homolog', 'composite'), 'divergence'),
        pch = c(24, 25, 24)
    )
    legend(
        0.54, 13.5, bty = 'n', cex = ifelse(presentation, 1.6, 1.2),
        rep(' ', 3),
        pch = c(NA, NA, 25)
    )
}

dev.off()














