#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

############
# Get data #
############

source('scripts/R_functions/load_genome_table.R')

presentation <- ifelse( "--presentation" %in% args, T, F)

genome_tab <- load_genome_table(c(1:5,15:18))
genome_tab <- genome_tab[nrow(genome_tab):1,]

# kick out those that have no TEs estimated
genome_tab <- genome_tab[!is.na(genome_tab$TEs), ]

# TO COMPARE WITH REPORTED TEs
# library('gsheet')
# literature_data <- read.csv(text = gsheet2text("https://docs.google.com/spreadsheets/d/1T4BHQxzMGMlWiNJ7G9OLPXzdNBCqMw0mSO7C_ggGEbc/edit?usp=sharing", format='csv'),
#                             stringsAsFactors = F, skip = 1, header = T, check.names = F)
# reported_TEs <- literature_data[,c("code","Repeats [ % ]", "TEs [ % ]")]
# colnames(reported_TEs) <- c('code', 'reported_rep', 'reported_TEs')
# not_analyzed <- reported_TEs[!(reported_TEs$code %in% rownames(genome_tab)),]
# reported_TEs <- reported_TEs[(reported_TEs$code %in% rownames(genome_tab)),]
# desired_order <- genome_tab$code
# genome_tab <- merge(genome_tab,reported_TEs)
#
# row.names(genome_tab) <- genome_tab$code
# genome_tab <- genome_tab[desired_order, 5:10]

source('scripts/R_functions/make_data_frame.R')

load_TE_tab <- function(x){
    file <- paste0('data/', x, '/dnaPipeTE/Counts.txt')
    if ( file.exists(file) ){
        file <- read.table(file)
        TE_table <- make_data_frame(as.character(file$V1[1:6]))
        TE_table[1,] <- file$V2[1:6] / file$V2[15]
        return(TE_table)
    }
    rep(NA, 6)
}

TEs <- as.data.frame(t(sapply(row.names(genome_tab), load_TE_tab))[,c(4,2,1,3,6)])
colnames(TEs) <- c('DNA', 'LINE', 'LTR', 'SINE', 'Helitron')
TEs <- TEs[,5:1] * 100
rownames(TEs)

# c('DNA', 'LINE', 'LTR', 'SINE', 'Helitron')
if( presentation ){
    TEs['5_Tge',] <- rev(c(11.57, 4.57, 3.05, 0.03, 0.39))
    TEs['4_Tte',] <- rev(c(11.65, 5.68, 3.27, 0.13, 0.54))
    TEs['3_Tms',] <- rev(c(12.67, 5.73, 3.14, 1.06, 0.45))
    TEs['2_Tsi',] <- rev(c(12.5, 4.57, 3.67, 1.07, 0.36))
    TEs['1_Tdi',] <- rev(c(10.72, 5.32, 3.08, 0.91, 0.46))
    TEs <- TEs[c(1:17, 30:34, 18:29),]
}


########
# plot #
########

pal <- c('#FF6A42', # DNA
         '#251792', # LINE
         '#65BD61', # LTR
         '#B966F4', # SINE
         '#8b0000') # helitron

if ( presentation ){
    pal <- c(grey(0.5),
             grey(0.2),
             grey(0.8),
             grey(0.35),
             grey(0.65))
}

gap_size <- ifelse(presentation, 1.5, 1)
spacing <- ifelse(presentation, 0.2, 0.15)
arthropods <- ifelse(presentation, 15, 10)
spaces <- c(rep(spacing, 2), gap_size,
            rep(spacing, 10), gap_size,
            rep(spacing, arthropods), gap_size,
            rep(spacing, 3), gap_size)
plot_bars <- function(index){
    if (index == 5){
        bar_sizes <- TEs[,5]
    } else {
        bar_sizes <- rowSums(TEs[,index:5])
    }
    barplot(bar_sizes, space = spaces,
            col = pal[6 - index],
            horiz = T, add = T, axes = F, names.arg = F)
}

bp <- barplot(rowSums(TEs), space = spaces,
              horiz = T, plot = F)
ylim <- range(bp)
xlim <- c(0,max(rowSums(TEs), na.rm = T))

# labels
names_split <- strsplit(genome_tab$species, "_")
genus_names <- sapply(names_split, function(x){substr(x[1], 1, 1)})
species_names <- sapply(names_split, function(x){ x[2] })
sp_labels <- paste(genus_names, species_names, sep = '. ')
sp_labels <- lapply(sp_labels, function(x){bquote(italic(.(x)))})

filename <- paste0('figures/',
                   ifelse(presentation, 'presentation/', ''),
                   'fig4_TEs',
                   ifelse(presentation, '_presentation', '') ,
                   '.pdf')
pdf(filename, width = ifelse(presentation, 9, 12), height = 8)

par(mar = c(5, 7, 1, 0.5) + 0.1)

plot(NULL, bty = 'n', axes=FALSE, xlim = xlim, ylim = ylim,
    xlab = "TEs [ % ]", ylab = NA, cex.lab = ifelse(presentation, 1.6, 1.3),
)
axis(1, cex.axis = ifelse(presentation, 1.4, 1.2))
if ( !presentation ){
    mtext(do.call(expression, sp_labels), line = 2.8, side = 2, at = bp, las = 2, adj = 0.5)
}
# axis(2, at = bp, labels, las = 1)

sapply(1:5, plot_bars)
if ( presentation ){
    legend('bottomright', bty = 'n', c('DNA', 'LINE', 'LTR', 'SINE', 'Helitron'),
           pch = 15, col = pal, cex = 1.4, horiz = T)
} else {
    legend('topright', bty = 'n', c('DNA', 'LINE', 'LTR', 'SINE', 'Helitron'),
           pch = 20, col = pal, cex = 1.2, title = 'TE classes')
}

dev.off()