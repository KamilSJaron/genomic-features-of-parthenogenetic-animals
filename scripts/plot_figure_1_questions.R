#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

require('gsheet')
require('RColorBrewer')

source('scripts/R_functions/wrap_labels.R')

presentation <- ifelse( "--presentation" %in% args, T, F)
refs <- ifelse( "--refs" %in% args, T, F)
both <- ifelse( "--both" %in% args, T, F)
tricolor <- ifelse( "--tricolor" %in% args, T, F)

pal <- c("white",brewer.pal(3,'BuGn')[-1], "grey") #[c(2,3,1)]#[c(2,3,5)]
shift <- 15

#################
# Proccess data #
#################

### BACKGROUND MATRIX ###

question_tab <- read.csv(text=gsheet2text("https://docs.google.com/spreadsheets/export?id=1T4BHQxzMGMlWiNJ7G9OLPXzdNBCqMw0mSO7C_ggGEbc&format=csv&gid=1327405223", format='csv'), stringsAsFactors=FALSE)

plotted_cols <- c(2:ncol(question_tab))
question_matrix <- as.matrix(question_tab[,plotted_cols])

convertor <- function(x){
  x <- substr(x, 1, 1)
    if(x == "b"){
        return(3)
    } else if(x == "q"){
        return(2)
    } else if(x == "d"){
    return(4)
  } else {
        return(0)
    }
}

# sapply(question_tab$species, function(x){expression(italic(x))})

# create a matrix of colours
heat_matrix <- matrix(sapply(question_matrix, FUN = convertor), nrow = nrow(question_tab))
# transpose and reverse for image
heat_matrix <- apply(heat_matrix, 2, rev)
heat_matrix <- t(heat_matrix)

if ( tricolor ){
  heat_matrix[heat_matrix == 3] <- 2
}

if ( refs ){
  ### REFERENCES ###
  ref_matrix <- matrix(sapply(question_matrix, FUN = function(x){ substr(x, 2, max(2, nchar(x))) }), nrow = nrow(question_tab))
  # separate references by , and take unique list
  ref_list <- unique(unlist(strsplit(t(ref_matrix), ',')))
  ref_tags <- lapply(1:length(ref_list), function(x){ x + shift})
  names(ref_tags) <- ref_list

  # substitute ref with a key
  for(ref in ref_list){
    ref_matrix <- gsub(ref, ref_tags[ref], ref_matrix)
  }

  # nicer formating ( [] around and space after comma)
  ref_matrix <- matrix(sapply(ref_matrix, FUN = function(x){ if(x != ""){ paste0('[', x, ']') } else { return("") } }), nrow = nrow(ref_matrix))
  ref_matrix <- gsub(",", ", ", ref_matrix)
}
if ( both | !refs ) {
  ### VALUES ###
  columns <- c('code',
               'increased mutation accomulation', # mutation accumulation
               'adaptive evolution',
               'haplotype_divergence[%]',
               'palindromes [ # ]',
               'gene conversion [ events / (generation * site) ]',
               'TEs [ % ]',
               'HGT [ % ]',
               'expansion of gene families',
               'missing sex related genes')

  literature_numbers <- read.csv(text = gsheet2text("https://docs.google.com/spreadsheets/d/1T4BHQxzMGMlWiNJ7G9OLPXzdNBCqMw0mSO7C_ggGEbc/edit?usp=sharing", format='csv'),
                                 stringsAsFactors = F, skip = 1, header = T, check.names = F)
  literature_numbers <- literature_numbers[, columns]

  to_round <- "Dcor1" != literature_numbers$code
  literature_numbers[to_round, "TEs [ % ]"] <- round(literature_numbers[to_round, "TEs [ % ]"], 1)

  sp_code <- function(code){ substr(code, 1, 4) }
  ### some of the species have more values -> range or list??
  squashed_lit_nums <- data.frame(code = unique(sp_code(literature_numbers$code)), stringsAsFactors = F)
  squashed_lit_nums[, columns[-1]] <- NA

  # remove "Ps59", "Ps79"
  squashed_lit_nums <- squashed_lit_nums[!squashed_lit_nums$code %in% c("Ps59", "Ps79"),]
  row.names(squashed_lit_nums) <- squashed_lit_nums$code

  merge_vals <- function(one_col_vals){
    # remove missing vals
    one_col_vals <- one_col_vals[!(is.na(one_col_vals) | one_col_vals == 'NA')]
    if ( length(one_col_vals) == 0 ){
      return("")
    }
    if ( length(one_col_vals) == 1 ){
      return(one_col_vals)
    }
    paste(one_col_vals, collapse = ', ')
  }

  for ( line in 1:nrow(squashed_lit_nums) ){
    sp <- sp_code(squashed_lit_nums[line, 'code'])
    if ( sp == "Pdav" ){
      to_squash <- sp_code(literature_numbers$code) %in% c("Pdav", "Ps59", "Ps79")
    } else {
      to_squash <- sp_code(literature_numbers$code) == sp
    }
    vals <- literature_numbers[to_squash, columns]
    if ( nrow(vals) > 1 ){
      squashed_lit_nums[sp, columns[-1]] <- apply(vals[,-1], 2, merge_vals)
    } else {
      vals[is.na(vals)] <- ""
      squashed_lit_nums[sp, columns[-1]] <- vals[-1]
    }
  }
}

########
# plot #
########

file_to_save <- paste0('figures/fig1_genomic_studies',
                       ifelse(refs, '_refs', ''),
                       ifelse(both, '_both', ''),
                       ifelse(presentation, '_presentation', ''), '.pdf')
# , family="Arial"
# height could be (ncol(ref_matrix) * 0.15 + 0.4) (which matches approximatelly the spacing required for new column)
if ( both ){
  height <- (ncol(ref_matrix) * 0.8 + 0.4)
} else {
  height <- 4.3
}

if (refs && both && tricolor){
    file_to_save <- "figures/SM_Figure_3_genomic_studies.pdf"
}

pdf(file_to_save, width = 8, height = height, pointsize = 8)

# I will need space on the left side of the plot and above
par(mar = c(0, 8, 4, 0) + 0.1)

# this will plot the hearmap but also keep the scale of the image (columns 1 .. questions and rows 1 .. sequenced genomes)
image(seq(dim(heat_matrix)[1]), seq(dim(heat_matrix)[2]), heat_matrix, col=pal, xaxt="n", yaxt="n", ylab="", xlab="")

# species names, the space in the end of each sp names is just for nice alignment
# These are not embeded into Illustartor picture but they can be regenerated by uncommenting following two lines
sp_labels <- lapply(paste0(question_tab$species, " "), function(x){bquote(italic(.(x)))})
mtext(do.call(expression, sp_labels), side = 2, at = rev(1:nrow(question_tab)), las = 2, adj = 0.6, line = 3.25)
# mtext(do.call(expression, sp_labels), side = 2, at = bp)

# names are taken from columns, but with spaces instead of periods
topics <- gsub("\\.", " ", colnames(question_tab)[plotted_cols])
if (presentation) {
  text(1:(ncol(question_tab) - 1), par("usr")[4] + 1.55 + c(0, 1.6), wrap.labels(topics, 10), xpd = TRUE)
} else{
  topics <- paste(topics, c('','', ' [%]', '[#]', '[e/(g*nt)]', '[%]', '[%]', '', ' [l/s] '))
  if ( both ){
    text(1:(ncol(question_tab) - 1), par("usr")[4] + 1.05, wrap.labels(topics, 11), xpd = TRUE)
  } else {
    text(1:(ncol(question_tab) - 1), par("usr")[4] + 2.1, wrap.labels(topics, 11), xpd = TRUE)
    # text(1:(ncol(question_tab) - 1) - 0.12, par("usr")[4] + 2.1, wrap.labels(topics, 10), xpd = TRUE, srt = 27)
  }
}

# create a black box around
box()
# lines(c(3.5, 3.5), c(-1000,1000), lwd = 1.5, lty = 1)
# lines(c(10.5, 10.5), c(-1000,1000), lwd = 1.5, lty = 2)

# adding references -> replace with numbers!!!
if (!presentation){
  if ( refs & !both ){
    for(line in nrow(ref_matrix):1){
      text(1:ncol(ref_matrix), line - 0.7, ref_matrix[nrow(ref_matrix) + 1 - line,], pos = 3, cex = 0.875)
    }
  }
  if ( !both | !refs ) {
    squashed_lit_nums <- squashed_lit_nums[,-1] # remove column with sp names
    for(line in nrow(squashed_lit_nums):1){
      text(1:ncol(squashed_lit_nums), line - 0.7, squashed_lit_nums[nrow(squashed_lit_nums) + 1 - line,], pos = 3, cex = 0.875)
    }
  }
  if ( both ){
    squashed_lit_nums <- squashed_lit_nums[,-1] # remove column with sp names
    for(line in nrow(squashed_lit_nums):1){
      text(1:ncol(squashed_lit_nums), line - 0.1, squashed_lit_nums[nrow(squashed_lit_nums) + 1 - line,], pos = 3, cex = 0.875)
      text(1:ncol(ref_matrix), line - 0.51, ref_matrix[nrow(ref_matrix) + 1 - line,], pos = 3, cex = 0.875)
    }
  }
}

dev.off()