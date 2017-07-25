#!/usr/bin/Rscript

# Enrichment of Mutations in and around Motifs 
# Tahmid Mehdi, Princess Margaret Cancer Centre - UHN, Dec 21, 2015

args <- commandArgs(trailingOnly = TRUE)
workingDir <- args[1]
mutations <- args[2] # bed file with mutation locations
motifs <- args[3] # bed file with motif locations
left <- as.numeric(args[4]) # length of left flanking sequence
right <- as.numeric(args[5]) # length of right flanking sequence
imgName <- args[6] # name of output image
imgWidth <- as.numeric(args[7]) # image width
imgHeight <- as.numeric(args[8]) # image height
tf <- args[9] # transcription factor whose motif is being plotted
binWidth <- as.numeric(args[10]) # width of the bins

setwd(workingDir)
library("GenomicRanges")
library("ggplot2")
library("stringr")

# From http://davetang.org/muse/2015/02/04/bed-granges/
bed_to_granges <- function(file){
  df <- read.table(file,
                   header=F,
                   stringsAsFactors=F)[,1:4]
  
  if(length(df) > 6){
    df <- df[,-c(7:length(df))]
  }
  
  if(length(df)<3){
    stop("File has less than 3 columns")
  }
  
  header <- c('chr','start','end','strand','id','score')
  names(df) <- header[1:length(names(df))]
  
  if('strand' %in% colnames(df)){
    df$strand <- gsub(pattern="[^+-]+", replacement = '*', x = df$strand)
  }
  
  if(length(df)==3){
    gr <- with(df, GRanges(chr, IRanges(start, end)))
  } else if (length(df)==4){
    gr <- with(df, GRanges(chr, IRanges(start, end), strand=strand))
  } else if (length(df)==5){
    gr <- with(df, GRanges(chr, IRanges(start, end), id=id, score=score))
  } else if (length(df)==6){
    gr <- with(df, GRanges(chr, IRanges(start, end), id=id, score=score, strand=strand))
  }
  return(gr)
}
# like bed_to_granges but for mutation files
mut_to_granges <- function(file){
  df <- read.table(file,
                   header=F,
                   stringsAsFactors=F)[,1:4]
  
  if(length(df) > 6){
    df <- df[,-c(7:length(df))]
  }
  
  if(length(df)<3){
    stop("File has less than 3 columns")
  }
  
  header <- c('chr','start','end','id','strand','score')
  names(df) <- header[1:length(names(df))]
  
  if('strand' %in% colnames(df)){
    df$strand <- gsub(pattern="[^+-]+", replacement = '*', x = df$strand)
  }
  
  if(length(df)==3){
    gr <- with(df, GRanges(chr, IRanges(start, end)))
  } else if (length(df)==4){
    gr <- with(df, GRanges(chr, IRanges(start, end), id=id))
  } else if (length(df)==5){
    gr <- with(df, GRanges(chr, IRanges(start, end), id=id, score=score))
  } else if (length(df)==6){
    gr <- with(df, GRanges(chr, IRanges(start, end), id=id, score=score, strand=strand))
  }
  return(gr)
}
# converts substitution types to uppercase & reverse complements if neccessary
convertSubType = function(type) {
  type <- toupper(type)
  if (type=="A>C") {return("T>G")
  } else if (type=="A>G") {return("T>C")
  } else if (type=="A>T") {return("T>A")
  } else if (type=="G>A") {return("C>T")
  } else if (type=="G>C") {return("C>G")
  } else if (type=="G>T") {return("C>A")
  } else {return(type)}
}
convertSubTypes = Vectorize(convertSubType)
# create GRanges objects
mut <- mut_to_granges(mutations)
mot <- bed_to_granges(motifs)
motWidth <- width(mot)[1]
# extend motifs 
flankedMot <- punion(punion(mot, flank(mot,left)), flank(mot,right,start=FALSE))
flankedMotWidth <- width(flankedMot)[1]
endMot <- flankedMotWidth - right + 0.5
mut <- mut[unique(queryHits(findOverlaps(mut, flankedMot)))]
# get nearest motif neighbours for each mutation
nn <- nearest(mut,mot,ignore.strand=TRUE)
# calculate distances of mutations from 5' ends of motif flanking regions & track substitution types
mutInMot <- findOverlaps(mut, mot)
if (length(mutInMot)==0) {
  distances <- c()
  types <- c()
} else {
  distances <- sapply(1:length(mutInMot), 
                      function(i) if(as.character(strand(flankedMot)[subjectHits(mutInMot)[i]])=="+") {
                        start(mut)[queryHits(mutInMot)[i]] - start(flankedMot)[subjectHits(mutInMot)[i]] + 1
                      } else { end(flankedMot)[subjectHits(mutInMot)[i]] - start(mut)[queryHits(mutInMot)[i]] + 1 })
  types <- mut[queryHits(mutInMot)]$id
}

mutNotInMot <- which(!(1:length(mut) %in% queryHits(mutInMot)))
distances <- c(distances, sapply(mutNotInMot, 
                                 function(i) if(as.character(strand(flankedMot)[nn[i]])=="+") {
                                   start(mut)[i] - start(flankedMot)[nn[i]] + 1
                                 } else { end(flankedMot)[nn[i]] - start(mut)[i] + 1 }))
types <- c(types, mut[mutNotInMot]$id)
types <- unname(convertSubTypes(types))
# create data frame with distances and types
mutDist <- data.frame(distances, types)
# filter distances
mutDist <- mutDist[1<=mutDist$distances & mutDist$distances<=flankedMotWidth, ]
# colours for histograms & streamgraphs
z <- c("deepskyblue2", "#252525", "firebrick2", "gray", "chartreuse3", "rosybrown1")
# sum up mutations for each distance-type combination
mutDist$freq <- 1
mutationDist <- aggregate(freq ~ distances + types, data = mutDist, sum)

#plot.stacked makes a stacked plot where each y series is plotted on top
#of the each other using filled polygons
#
#Arguments include:
#'x' - a vector of values
#'y' - a matrix of data series (columns) corresponding to x
#'order.method' = c("as.is", "max", "first") 
#  "as.is" - plot in order of y column
#  "max" - plot in order of when each y series reaches maximum value
#  "first" - plot in order of when each y series first value > 0
#'col' - fill colors for polygons corresponding to y columns (will recycle)
#'border' - border colors for polygons corresponding to y columns (will recycle) (see ?polygon for details)
#'lwd' - border line width for polygons corresponding to y columns (will recycle)
#'...' - other plot arguments

plot.stacked <- function(
  x, y, 
  order.method = "as.is",
  ylab="", xlab="", 
  border = NULL, lwd=1, 
  col=rainbow(length(y[1,])),
  ylim=NULL, bp, left, right,
  ...
){
  
  if(sum(y < 0) > 0) error("y cannot contain negative numbers")
  
  if(is.null(border)) border <- par("fg")
  border <- as.vector(matrix(border, nrow=ncol(y), ncol=1))
  col <- as.vector(matrix(col, nrow=ncol(y), ncol=1))
  lwd <- as.vector(matrix(lwd, nrow=ncol(y), ncol=1))
  
  if(order.method == "max") {
    ord <- order(apply(y, 2, which.max))
    y <- y[, ord]
    col <- col[ord]
    border <- border[ord]
  }
  
  if(order.method == "first") {
    ord <- order(apply(y, 2, function(x) min(which(x>0))))
    y <- y[, ord]
    col <- col[ord]
    border <- border[ord]
  }
  
  top.old <- x*0
  polys <- vector(mode="list", ncol(y))
  for(i in seq(polys)){
    top.new <- top.old + y[,i]
    polys[[i]] <- list(x=c(x, rev(x)), y=c(top.old, rev(top.new)))
    top.old <- top.new
  }
  
  if(is.null(ylim)) ylim <- range(sapply(polys, function(x) range(x$y, na.rm=TRUE)), na.rm=TRUE)
  plot(x,y[,1], ylab=ylab, xlab=xlab, ylim=ylim, t="n", xaxt="n", bty="n", ...)
  axis(1, at=c(1,bp), labels=c(-left, right))
  for(i in seq(polys)){
    polygon(polys[[i]], border=border[i], col=col[i], lwd=lwd[i])
  }
  
}
# create a matrix where each row corresponds to a bp & columns are types
M <- matrix(0, nrow=flankedMotWidth, ncol=6)
colnames(M) <- c("C>A","C>G","C>T","T>A","T>C","T>G")
# the ith row, jth column is the number of mutations in the ith bp & jth type
for (i in 1:nrow(mutationDist)) {
  M[mutationDist$distances[i], mutationDist$types[i]] <- mutationDist$freq[i]
}
# aggregate the rows in matrix based on binwidth
B <- matrix(0, nrow=0, ncol=6)
if (binWidth==1) {
  B <- M
} else { # sum number of mutations for each type for a bin
  for (i in seq(1,flankedMotWidth,binWidth)) {
    if (i==flankedMotWidth) {
      B <- rbind(B, M[i,])
    } else if (i+binWidth-1 <= flankedMotWidth) {
      B <- rbind(B, colSums(M[i:(i+binWidth-1),]))
    } else {
      B <- rbind(B, colSums(M[i:flankedMotWidth,]))
    }
  }
}
# generate plot
tiff(paste(imgName,"_stream.tif", sep=""), width=imgWidth, height=imgHeight, res=100)
layout(matrix(1:2, ncol=2), width=c(6,1), height=3.25, respect=TRUE)
par(mar=c(2,5,0,1), cex=0.75)
plot.stacked(seq(nrow(B)),B, yaxs="i", col=z, border="white", lwd=0.5, bp=nrow(B), left=left, right=right)
abline(v=left/binWidth, lty=2, lwd=2)
mtext(paste("N=",length(mot)," motif occurences", sep=""), line=1, side=3)
mtext(paste(tf," + flanking region (bp)", sep=""), line=2.5, side=1)
mtext("Number of Mutations", line=2.5, side=2)
# legend
par(mar=c(1,1,1,1), cex=0.75)
plot(1,t="n", xlab="", ylab="", axes=FALSE)
legend("left", legend=c("C>A","C>G","C>T","T>A","T>C","T>G"),
       border="white", fill=z, bty="n", cex=2, y.intersp=2)
dev.off()


# generate plot
pdf(paste(imgName,"_stream.pdf", sep=""), width=4, height=3, pointsize = 14 )
layout(matrix(1:2, ncol=2), width=c(6,1), height=3.25, respect=TRUE)
par(mar=c(2,5,0,1), cex=0.75)
plot.stacked(seq(nrow(B)),B, yaxs="i", col=z, border="white", lwd=0.5, bp=nrow(B), left=left, right=right)
abline(v=left/binWidth, lty=2, lwd=2)
abline(v=(left+motWidth+binWidth)/binWidth, lty=2,lwd=2)
mtext(paste("N=",length(mot)," motif occurences", sep=""), line=1, side=3)
mtext(paste(tf," + flanking region (bp)", sep=""), line=2.5, side=1)
mtext("Number of Mutations", line=2.5, side=2)
# legend
par(mar=c(1,1,1,1), cex=0.75)
plot(1,t="n", xlab="", ylab="", axes=FALSE)
legend("left", legend=c("C>A","C>G","C>T","T>A","T>C","T>G"),
       border="white", fill=z, bty="n", cex=2, y.intersp=2)
dev.off()


#for (subType in colnames(M)) {
#  fit <- smooth.spline(M[,subType] ~ seq(nrow(M)), df=20)
#  M[,subType] <- fit$y
#}
#M[M<0] <- 0

