library('shape')
library('grid')

# empty cavas
emptyplot <- function(){
    plot(NULL, bty = 'n', axes=FALSE, xlim = c(0, 1), ylim = c(0,1), xlab = '', ylab = '')
}

getCenter <- function(x, ch){
    y_center <- 0.5
    x_center <- ((x - 1) / 3) + sm + (p / 30) * ((8 * ch) - 1)
    c(x_center, y_center)
}

plotMutations <- function(x, .mutation_matrix){
    # then draw grey boxes bellow heterozygous positions
    # then draw mutations
    num_of_ch <- ncol(.mutation_matrix)
    heterozygous_loci <- !apply(mutation_matrix, 1, function(x){ all(x == x[1])})

    locus_base <- ((getCenter(x, 1) + getCenter(x, num_of_ch)) / 2) - c(0, (ch_height / 2))


    for(pos in 1:loci){
        if(heterozygous_loci[pos]){
            grid.rect(locus_base[1],
                      locus_base[2] + pos * (ch_height / loci),
                      width = (p/5 * num_of_ch) + (p/15 * (num_of_ch - 1)),
                      height = ch_height / loci,
                      gp=gpar(col=NA, fill="lightgrey"))
        }
    }

    for(.ch in 1:num_of_ch){
        ch_base <- getCenter(x,.ch) - c(0, ch_height / 2)
        for(pos in 1:loci){
            .mut = .mutation_matrix[pos,.ch]
            if(.mut != ""){
                grid.rect(ch_base[1],
                          ch_base[2] + pos * (ch_height / loci),
                          width = p / 5,
                          height = ch_height / loci,
                          gp=gpar(col=NA, fill=.mut))
            }
        }
    }

    mean(heterozygous_loci)
}

mut_pal <- c('red', 'green', 'blue', 'yellow')

mutate <- function(locus_vector, f_mut){
    mutations <- rbinom(length(locus_vector), 1, f_mut)
    locus_vector[which(mutations == 1)] <- sample(mut_pal, sum(mutations), replace = T)
    locus_vector
}

sm <- 0.03 # an extra margin on right within a box
p <- 0.22 # proportion of a box
ch_height <- 4/5 # chromosome height
loci <- 100

mut_matrices <- list()
mut_matrices[[1]] <- sapply(rep(0.02, 2), function(x){ mutate(rep("", loci), x)})
mut_matrices[[2]] <- sapply(rep(0.015, 3), function(x){ mutate(rep("", loci), x)})

A <- mutate(rep("", loci), 0.01)
A1 <- mutate(A, 0.005)
A2 <- mutate(A, 0.005)

B <- mutate(rep("", loci), 0.01)
B1 <- mutate(B, 0.005)
B2 <- mutate(B, 0.005)

mut_matrices[[3]] <- matrix(c(A1, A2, B1, B2), ncol = 4)

pdf('figures/box3_scheme.pdf')
emptyplot()

for(x in c(1,2,3)){
    mutation_matrix <- mut_matrices[[x]]
    heterozygosity <- plotMutations(x, mutation_matrix)
    text_pos <- getCenter(x, 1 + (0.5 * x))
    # grid.text(round(heterozygosity, 2), text_pos[1], 0.03)
    for(ch in 1:(x+1)){
        ch_center <- getCenter(x,ch)
        grid.roundrect(ch_center[1], ch_center[2], width = p/5, height = ch_height, gp=gpar(fill="#00000000"))
    }
}

grid.text("diploid", getCenter(1,1.5)[1], 0.05, gp=gpar(fontsize=20))
grid.text("triploid", getCenter(2,2)[1], 0.05, gp=gpar(fontsize=20))
grid.text("tetraploid", getCenter(3,2.5)[1], 0.05, gp=gpar(fontsize=20))

dev.off()
