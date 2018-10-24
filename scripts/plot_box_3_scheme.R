library('shape')
library('grid')

# empty cavas
emptyplot <- function(){
    plot(NULL, bty = 'n', axes=FALSE, xlim = c(0, 1), ylim = c(0,1), xlab = '', ylab = '')
}

getCenter <- function(x, y, ch){
    y_center <- (y / 4) - ch_height / 2 - 0.01
    x_center <- ((x - 1) / 2) + sm + (p / 30) * ((8 * ch) - 1) - 0.15 * (x - 1)
    c(x_center, y_center)
}

plotMutations <- function(x, y, .mutation_matrix){
    # then draw grey boxes bellow heterozygous positions
    # then draw mutations
    num_of_ch <- ncol(.mutation_matrix)
    heterozygous_loci <- !apply(mutation_matrix, 1, function(x){ all(x == x[1])})

    locus_base <- ((getCenter(x,y,1) + getCenter(x,y,num_of_ch)) / 2) - c(0, (ch_height / 2))


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
        ch_base <- getCenter(x,y,.ch) - c(0, ch_height / 2)
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

sm <- 0.05 # an extra margin on right within a box
p <- 0.15 # proportion of a box
ch_height <- 1/5 # chromosome height
loci <- 50

mut_matrices <- lapply(1:2,function(x){list()})
mut_matrices[[1]][[1]] <- matrix(mutate(rep("", loci), 0), ncol = 1)
mut_matrices[[2]][[1]] <- matrix(mutate(rep("", loci), 0.05), ncol = 1)
mut_matrices[[1]][[2]] <- sapply(rep(0.02, 2), function(x){ mutate(rep("", loci), x)})
mut_matrices[[2]][[2]] <- sapply(rep(0.08, 2), function(x){ mutate(rep("", loci), x)})
mut_matrices[[1]][[3]] <- sapply(rep(0.02, 3), function(x){ mutate(rep("", loci), x)})
mut_matrices[[2]][[3]] <- sapply(rep(0.08, 3), function(x){ mutate(rep("", loci), x)})

A <- mutate(rep("", loci), 0.03)
A1 <- mutate(A, 0.01)
A2 <- mutate(A, 0.01)

B <- mutate(rep("", loci), 0.03)
B1 <- mutate(B, 0.01)
B2 <- mutate(B, 0.01)

mut_matrices[[1]][[4]] <- matrix(c(A1, A2, B1, B2), ncol = 4)

A <- mutate(rep("", loci), 0.03)
A1 <- mutate(A, 0.01)
A2 <- mutate(A, 0.01)
A3 <- mutate(A, 0.01)
B <- mutate(rep("", loci), 0.05)
mut_matrices[[2]][[4]] <- matrix(c(A1, A2, A3, B), ncol = 4)

pdf('figures/box3_scheme.pdf')
emptyplot()

for(x in c(1,2)){
    for(y in c(1,2,3,4)){
        mutation_matrix <- mut_matrices[[x]][[5 - y]]
        heterozygosity <- plotMutations(x, y, mutation_matrix)
        text_pos <- getCenter(x,y,6)
        if( y != 4){
            grid.text(round(heterozygosity, 2), text_pos[1], text_pos[2])
        }
        for(ch in 1:(5-y)){
            ch_center <- getCenter(x,y,ch)
            grid.roundrect(ch_center[1], ch_center[2], width = p/5, height = ch_height, gp=gpar(fill="#00000000"))
        }
    }
}

grid.text("the ancestral chromosome", 0.285, getCenter(1,4,1)[2])
grid.text("a chromosome with mutations", 0.635, getCenter(1,4,1)[2])

grid.text("heterozygosity is the proportion of sites with more than one allele", 0.4, getCenter(1,4,1)[2] - 0.125)

grid.text("tetraploid", 0.03, getCenter(1,1,1)[2], rot = 90)
grid.text("triploid", 0.03, getCenter(1,2,1)[2], rot = 90)
grid.text("diploid", 0.03, getCenter(1,3,1)[2], rot = 90)

dev.off()
