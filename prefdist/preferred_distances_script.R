#!/usr/bin/Rscript

# INPUT args 
# peak fasta file
# chr1 fasta file for background
# anchor motif
# other motif

# OUTPUT args
# table_file
# image_pdf

args<-commandArgs(TRUE)
library("Biostrings")
library("stringr")
source("preferred_distances4_source.R")
#####################################################
# Get the data ready for preferred distance cmputation
#####################################################
#mot<-read.delim(file="tal1_ecfc_motifs.csv")
# Extract  motif sequences
#motifs<-rownames(mot);
#input
#peak_fasta_file<-args[1]
peak_fasta_file="tal1_ecfc_peaks.fasta"

#chr1_fasta_file<-args[2]
chr1_fasta_file<-"background_chr1/chr1c.fa"
#anchor_motif<-args[3]
anchor_motif<-"CANNTG"
#other_motif<-args[4]
other_motif<-"WGATAA"
#output
#table_file<-args[5]
table_file<-"flat.csv"
#image_file<-args[6]
image_file<-"pd.pdf"

# OLD Background  sequences 
#fasta<-read.DNAStringSet("tal1_ecfc_flank_bg.fasta",use.names=TRUE)

####################################################
#
# anchor motif, with respect to  which we 
# explore whether there is a preferred distance
# between it and other motifs
# We look for preferred distances in a set of peak  
# sequences and corresponding background sequences 
####################################################

################### Initialization ######################
# EBOX 
ret_rez_all<-list();
di<-seq(from=-20,to=20,by=1);

################### Composites ##########################
pat1<-str_trim(str_replace(str_replace(anchor_motif,"^[N]*","  "),"[N]*$","  "),side="both");
pat2<-str_trim(str_replace(str_replace(other_motif,"^[N]*","  "),"[N]*$","  "),side="both");
pat_o<-generateComposites(dist=di,pat1=consensusTOregexp(pat1),pat2=consensusTOregexp(pat2) )
################# PEAKS

######################################################
# Peak sequences 
fasta<-read.DNAStringSet(peak_fasta_file,use.names=TRUE)
peak_fasta<-extract_strings(fasta)

dim(peak_fasta)
rm(fasta)

ret_o<-count_composites3(data=peak_fasta,pat=pat_o)
# pattern occurrence matrix
flat=ret_o[[1]]
write.table(flat, file=table_file, sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)


# distance plot
rm(peak_fasta)
pd=ret_o[[2]]

##########BACKGROUND
#This is very time and resources consuming operation

#########################  Background
fasta<-read.DNAStringSet(chr1_fasta_file,use.names=TRUE)
backg_fasta<-extract_strings(fasta)
dim(backg_fasta)
rm(fasta); # clear the memory 

ret_b<-count_composites3(data=backg_fasta,pat=pat_o)
# we save the result and process it later with the result for peaks
pdb=ret_b[[2]]
rm(ret_b)

################## PLOT
# when we have pd and pdb , we plot the occurrences
npdb=sum(as.numeric(paste(pdb$lengths)))
npd=sum(as.numeric(paste(pd$lengths)))

# peaks
xp=as.numeric(paste(pd$values))
# yp=as.numeric(paste(m[,2]))/npd;
yp=as.numeric(paste(pd$lengths))/npd;

# background
xb=as.numeric(paste(pdb$values))
yb=as.numeric(paste(pdb$lengths))/npdb;

# zero distance in the background
yb[21]=yb[21]/2;
 
# plot
#plotname=paste(c(anchor_motif,other_motif,"_peaks_backgr.pdf"),collapse="")
mainpeaks=paste(c(anchor_motif,"-",other_motif," in peaks"),collapse="")
mainbackg=paste(c(anchor_motif,"-",other_motif," in Chromosome 1"),collapse="")

pdf(image_file,width=12,height=7)
par(mfcol=c(1,2))
 limit<-c(0,max(yp+0.001))
plot(xp,yp,type='h', main=mainpeaks,ylab='relative frequency', xlab='distance (nt)',ylim=limit);
plot(xb,yb,type='h', main=mainbackg,ylab='relative frequency', xlab='distance (nt)',ylim=limit) 
dev.off();
#############################################################################
#  END of Computations for preferred distances
###############################################################################