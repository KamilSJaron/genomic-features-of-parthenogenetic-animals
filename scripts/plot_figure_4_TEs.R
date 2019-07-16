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

########
# plot #
########

pal <- c('#FF6A42', # DNA
         '#251792', # LINE
         '#65BD61', # LTR
         '#B966F4', # SINE
         '#8b0000') # helitron

spaces <- c(0.1, 0.1, 1, rep(0.1, 10),1 , rep(0.15, 7), 1, rep(0.15, 3), 1)
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
genus_names <- sapply(names_split, substr, 1, 1)[1,]
species_names <- sapply(names_split, function(x){ x[2] })
species_names[10] <- c("davidi")
sp_labels <- paste(genus_names, species_names, sep = '. ')
sp_labels[9] <- "Panagrolaimus sp."
sp_labels <- lapply(sp_labels, function(x){bquote(italic(.(x)))})

filename <- paste0('figures/fig4_TEs',
                   ifelse(presentation, '_presentation', '') ,
                   '.pdf')
pdf(filename, width = 12, height = 8)

par(mar = c(5, 7, 1, 0.5) + 0.1)

plot(NULL, bty = 'n', axes=FALSE, xlim = xlim, ylim = ylim,
    xlab = "TEs [ % ]", ylab = NA, cex.lab = ifelse(presentation, 2, 1.3),
)
axis(1, cex.axis = ifelse(presentation, 1.6, 1.2))
if ( !presentation ){
    mtext(do.call(expression, sp_labels), line = 2.8, side = 2, at = bp, las = 2, adj = 0.5)
}
# axis(2, at = bp, labels, las = 1)

sapply(1:5, plot_bars)
legend('topright', bty = 'n', c('DNA', 'LINE', 'LTR', 'SINE', 'Helitron'),
       pch = 20, col = pal, cex = ifelse(presentation, 2, 1.2), title = 'TE classes')

dev.off()