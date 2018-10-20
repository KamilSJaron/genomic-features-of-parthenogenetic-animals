############
# Get data #
############

genome_tab <- read.table('tables/genome_table.tsv',
                         header = T, stringsAsFactors = F, skip = 1, check.names = F)
rownames(genome_tab) <- genome_tab$code
# take only relevant columns and reverse the order of rows
genome_tab <- genome_tab[nrow(genome_tab):1,c(1,3:5,9:12,15:18)]
# kick out those that have no
genome_tab <- genome_tab[!(is.na(genome_tab$TEs) & is.na(genome_tab$repeats) & is.na(genome_tab$complete)),]

hyb_origins <- c("no", "unknown", "yes")
repr_modes  <- c("gamete_duplication", "terminal_fusion", "central_fusion", "unknown_automixis",
                 "unknown", "functional_apomixis")

# ordering
genome_tab$hybrid_origin <- ordered(genome_tab$hybrid_origin, levels=hyb_origins )
genome_tab$reproduction_mode <- ordered(genome_tab$reproduction_mode, levels=repr_modes )

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

make_data_frame <- function (variables) {
    asm_template <- matrix(ncol = length(variables))
    colnames(asm_template) <- variables
    asm_template
}

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
TEs <- TEs[,5:1]

########
# plot #
########

pal <- c('#FF6A42', # DNA
         '#251792', # LINE
         '#65BD61', # LTR
         '#B966F4', # SINE
         '#8b0000') # helitron

spaces <- c(0.1, 0.1, 1, rep(0.1, 15),1 , rep(0.15, 8), 1, rep(0.15, 3), 1)
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

pdf('figures/fig3_TEs_BUSCO.pdf', width = 12, height = 8)
par(mfrow=c(1,2))
par(mar = c(5, 4, 1, 0.5) + 0.1)

plot(NULL, bty = 'n', axes=FALSE, xlim = xlim, ylim = ylim,
    xlab = "TEs [ % ]", ylab = NA, cex.lab = 1.3
)
axis(1)
axis(2, at = bp, row.names(TEs), las = 1)

sapply(1:5, plot_bars)
legend('topright', bty = 'n', c('DNA', 'LINE', 'LTR', 'SINE', 'Helitron'),
       pch = 20, col = pal, cex = 1, title = 'TE classes')

busco_pal <- c('#fed976', '#fd8d3c','grey')
legend('bottomright', bty = 'n', c('single', 'duplicated', 'missing'),
        pch = 20, col = busco_pal, cex = 1, title = 'BUSCO')

# plot BUSCO
# 'c(bottom, left, top, right)'
par(mar = c(5, 1, 1, 2) + 0.1)

barplot(c(rep(100, nrow(genome_tab) - 13), NA, NA, rep(100, 11)),
        xlab = "BUSCO [ % ]", cex.lab = 1.3,
        col = busco_pal[3], horiz = T, axes = F, names.arg = F, space = spaces)
axis(1)

barplot(genome_tab$complete + genome_tab$fragmented, space = spaces,
        col = busco_pal[2], horiz = T, add = T, axes = F, names.arg = F)
barplot(genome_tab$complete + genome_tab$fragmented - genome_tab$duplicated, space = spaces,
        col = busco_pal[1], horiz = T, add = T, axes = F, names.arg = F)

dev.off()