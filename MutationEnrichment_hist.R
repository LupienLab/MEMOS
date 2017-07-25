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


setwd("/Users/lupienlab3/Documents/tahmid/PCa")
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

# save image to current directory
pdf(paste(imgName,"_hist.pdf", sep=""), width=imgWidth, height=imgHeight, res=3000)
ggplot(mutDist, aes(distances, fill=types)) + geom_bar() +
  scale_x_continuous(breaks=c(1,flankedMotWidth), labels=c(-left, right)) + theme_bw() +
  theme(panel.grid.minor=element_blank(),panel.grid.major=element_blank(),
        axis.title.x=element_text(hjust=((left+motWidth/2)/flankedMotWidth))) +
  geom_vline(xintercept=left+0.5, lty=2) + geom_vline(xintercept=endMot, lty=2) + guides(fill=guide_legend(title=NULL)) +
  #  geom_vline(xintercept=left+11+0.5) + geom_vline(xintercept=endMot-11) +
  #  stat_smooth(method="loess", formula=df$freq~as.numeric(df$pos), size=1, aes(group=1), se=FALSE, span=0.25) +
  labs(x=paste(tf," + flanking region (bp)", sep=""), title=paste("N=",length(mot)," motif occurences", sep=""), y="Number of Mutations") +
  scale_fill_manual(values=z)
dev.off()
