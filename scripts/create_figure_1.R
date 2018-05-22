library('gsheet')
library('RColorBrewer')

pal <- brewer.pal(3,'BuGn')#[c(2,3,1)]#[c(2,3,5)]

# functions to wrap text from https://stackoverflow.com/a/20241729/2962344
# Core wrapping function
wrap.it <- function(x, len)
{
  sapply(x, function(y) paste(strwrap(y, len),
                              collapse = "\n"),
         USE.NAMES = FALSE)
}


# Call this function with a list or vector
wrap.labels <- function(x, len)
{
  if (is.list(x))
  {
    lapply(x, wrap.it, len)
  } else {
    wrap.it(x, len)
  }
}

############
# Get data #
############

question_tab <- read.csv(text=gsheet2text("https://docs.google.com/spreadsheets/export?id=1T4BHQxzMGMlWiNJ7G9OLPXzdNBCqMw0mSO7C_ggGEbc&format=csv&gid=1327405223", format='csv'), stringsAsFactors=FALSE)
question_tab <- question_tab[c(26, 27:30, 1, 3, 2, 20:22, 25, 4:5, 6:8, 19, 9:18, 23, 24),]

question_matrix <- as.matrix(question_tab[,c(3:ncol(question_tab))])

convertor <- function(x){
	if(x == "yes"){
		return(1)
	} else if(x == "no"){
		return(0)
	} else {
		return(0.5)
	}
}

# sapply(question_tab$species, function(x){expression(italic(x))})

question_matrix <- matrix(sapply(question_matrix, FUN = convertor), nrow = nrow(question_tab))
# transpose and reverse for image
question_matrix <- apply(question_matrix, 2, rev)
question_matrix <- t(question_matrix)

########
# plot #
########

pdf('figures/fig1_genomic_studies.pdf', width=12, height=8)

# I will need space on the left side of the plot and above
par(mar = c(0, 7.5, 5.5, 0) + 0.1)

# this will plot the hearmap but also keep the scale of the image (columns 1 .. questions and rows 1 .. sequenced genomes)
image(seq(dim(question_matrix)[1]), seq(dim(question_matrix)[2]), question_matrix,
col=pal, xaxt="n", yaxt="n", ylab="", xlab="")

# species names, the space in the end of each sp names is just for nice alignment
sp_labels <- lapply(paste0(question_tab$species, " "), function(x){bquote(italic(.(x)))})
mtext(do.call(expression, sp_labels), side = 2, at = rev(1:nrow(question_tab)), las = 1)

# names are taken from columns, but with spaces instead of periods
topics <- gsub("\\.", " ", colnames(question_tab)[3:ncol(question_tab)])
text(1:ncol(question_tab), par("usr")[4] + 2.3, wrap.labels(topics, 10), srt = 45, xpd = TRUE)

# create a black box around
box()

dev.off()