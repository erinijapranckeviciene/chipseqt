#preferred distances
# input 
# Peaks fasta
# anchor motif
# other motifs
# background sequence


#############################################################################################
#  Computations for preferred distances
#############################################################################################
# Libraries
library("Biostrings")
library("stringr")
source("preferred_distances4_source.R")
#####################################################
# Get the data ready for preferred distance cmputation
#####################################################
# load DE NOVO Motifs
mot<-read.delim(file="tal1_ecfc_motifs.csv")
# Extract  motif sequences
motifs<-rownames(mot);

# OLD Background  sequences 
#fasta<-read.DNAStringSet("tal1_ecfc_flank_bg.fasta",use.names=TRUE)


#############################

####################################################
#
# EBOX is an anchor motif, with respect to  which we 
# explore whether there is a preferred distance
# between EBOX and other motifs, known and 
# identified by the de novo search
#
# We look for preferred distances in a set of peak  
# sequences and corresponding background sequences 
# that comprise left and right flanking regions of the 
# peak region. Resulting matrixes are collected in a 
# list, the maximum ylim for plot
# is identified , plots are generated using the uniform scale   
####################################################

################### Initialization ######################
# EBOX 
ebox<-c("CANNTG")
ret_rez_all<-list();
m_names<-c();
di<-seq(from=-20,to=20,by=1);
################### Computation ##########################
# Known motifs ETS and other  MOTIFS THAT WE FOUND ( WGATAA among them)
#
motifs_all=motifs;
#VAGGAAR; AGGAAA; VCAGGA
motifs_all<-c("VAGGAAR","AGGAAA","VCAGGA",motifs)

for( i in 1:1)
{
m_name<-str_trim(str_replace(str_replace(motifs_all[i],"^[N]*","  "),"[N]*$","  "),side="both");
pat1=ebox;
pat2=m_name;
pat_o<-generateComposites(dist=di,pat1=consensusTOregexp(pat1),pat2=consensusTOregexp(pat2) )

################# PEAKS
######################################################
# Peak sequences 
fasta<-read.DNAStringSet("tal1_ecfc_peaks.fasta",use.names=TRUE)
peak_fasta<-extract_strings(fasta)
dim(peak_fasta)
rm(fasta)
ret_o<-count_composites3(data=peak_fasta,pat=pat_o)
# pattern occurrence matrix
#matches=ret_o[[1]]
# distance plot
# release memory
#save the result
fname=paste(c(m_name,"_ret_o.Rdata"),collapse="");
save(ret_o,fname,file=fname)
rm(peak_fasta)
pd=ret_o[[2]]

##########BACKGROUND
#This is very time and resources consuming operation

#########################  Background
fasta<-read.DNAStringSet("background_chr1/chr1c.fa",use.names=TRUE)
backg_fasta<-extract_strings(fasta)
dim(backg_fasta)
rm(fasta); # clear the memory 
ret_b<-count_composites3(data=backg_fasta,pat=pat_o)
# we save the result and process it later with the result for peaks
# save the ret_b background
fname=paste(c(m_name,"_ret_b.Rdata"),collapse="");
save(ret_b,fname,file=fname)
pdb=ret_b[[2]]

################## PLOT
# when we have pd and pdb , we plot the occurrences
npdb=sum(as.numeric(paste(pdb$lengths)))
 npd=sum(as.numeric(paste(pd$lengths)))
 # peaks
 xp=m[,1]; xp[21]=-0.0001;
# xp=as.numeric(paste(pd$values))
 yp=as.numeric(paste(m[,2]))/npd;
# yp=as.numeric(paste(pd$lengths))/npd;
 # background
 #xb=as.numeric(paste(pdb$values))
 xb=xp;
 yb=as.numeric(paste(m[,3]))/npdb;
 
 #yb=as.numeric(paste(pdb$lengths))/npdb;
 # zero distance in the background
#yb[21]=yb[21]/2;
 
# plot
plotname=paste(c("EBOX-",m_name,"_peaks_backgr.pdf"),collapse="")
mainpeaks=paste(c("EBOX-",m_name," in TAL1 ECFC peaks"),collapse="")
mainbackg=paste(c("EBOX-",m_name," in Chromosome 1"),collapse="")

pdf(plotname,width=12,height=7)
#pdf("EBOX-WGATAA_peaks_backgr.pdf",width=12,height=7)
 par(mfcol=c(1,2))
 limit<-c(0,max(yp+0.001))
plot(xp,yp,type='h', main=mainpeaks,ylab='relative frequency', xlab='distance (nt)',ylim=limit);
plot(xb,yb,type='h', main=mainbackg,ylab='relative frequency', xlab='distance (nt)',ylim=limit) 
 dev.off();
}



#save.image("preferred_distances_comp.Rdata")
# Save data
#write.table(flat, file="preferred_distances_of_composites.csv", sep="\t",quote=FALSE,row.names=FALSE)
#############################################################################################
#  END of Computations for preferred distances
#############################################################################################





# functions
