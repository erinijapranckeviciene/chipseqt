#!/data/binaries/R-3.0.1//bin/Rscript
# args motif fasta_file table_file
args <- commandArgs(TRUE)
#
#library("Biostrings")
#library("stringr")
source("/data/projects/ep_links/galaxy/import_scripts/functions.R")
#

###############################################################################################
# OCCURRENCES OF THE MOTIF in PEAK FASTA FILE
###############################################################################################
#
#
#  Functions in this Program
#  consensusTOregexp<-function(consstr="WGATAAN")
#

motif<-args[1];
fasta_file<-args[2];
table_file<-args[3];

# TEST SRIPT IN R CONSOLE 
#motif<-"WGATAA"
#fasta_file<-"tal1_ecfc_peaks.fasta"
#table_file<-"occ.csv"



outfile<-file(table_file,open="wt");

# redirect otput
sink(stdout(),type="message");
#sink(report,type="output");

fa<-readDNAStringSet(fasta_file,use.names=TRUE)
fstr<-str_trim(str_replace(str_replace(motif,"^[N]*","  "),"[N]*$","  "),side="both");
rstr<-reverse(fstr)

# make the single string into consensus str
forwstr<-consensusTOregexp(fstr)
revstr<-consensusTOregexp(rstr)

print(forwstr)
print(revstr)

#WE HAVE to transform DNA StringSet into the 
peaks<-extract_strings(fa)

#FOR EACH peak find the patterns
#return the list containing patterns, class, position
#the row in the data frame will contain
# peak , comma-separated patterns, comma-separated classes of the patterns, 
#        comma-separated positions of the patterns, 
# for the reverse pattern start should interpreted as stop, and stop as start

# get the occurrences of the motif in peaks

d<-patterns_in_peaks(peaks,forwstr,revstr)


write.table(d,file=outfile,sep="\t",quote=FALSE, col.names=TRUE, row.names=FALSE )
###############################################################################################
quit(status=0,save="no",runLast=FALSE);