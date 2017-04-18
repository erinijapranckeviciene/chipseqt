# For any list, given by symbols you can extract EntrezIDs by synergizer 
http://llama.mshri.on.ca/synergizer/translate/

# This is GO analysis of the list of gennes by GO stats package
# The ToppGene web Tool does it as well, but there you have only one gene universe- whole genome. 
# 

#load the library
library("GOstats");
# The genome, you use should be installed
# source("http://bioconductor.org/biocLite.R")
# biocLite("org.Hs.eg.db")
# org.Hs.eg.db in TAL1 ECFC case, this package takes time to install

# This example is based on Tal1 ECFC microarray
# Universe consists of all expressed genes in TAL1 ECFC experiment, that are on microarray
# up_down are the genes, that change significantly on the array


#Read in the data
#Genes of the Universe - column of EntrezID identifiers 
uni<-(read.delim("universe.ids"));

#Read the table with microarray result , fields separated by TAB
# In the table there already are entrez IDs 
# If the list of genes is given just by a columns of gene symbols, 
# You need to transform it to EntrezIDs, I do it with synergizer
# But there are means to do this in R as well. 
data<-read.delim("up_down.csv")


#########################################################
# Function you have to load into the workspace
#
#
##########################################################
GOenrichment<-function(ont, universeID, geneID)
##########################################################
{
	#ont<-"CC";
params<-new("GOHyperGParams",geneIds=geneID,universeGeneIds=universeID,annotation="org.Hs.eg",ontology=ont,pvalueCutoff=0.05,conditional=FALSE,testDirection="over");

hgOver<-hyperGTest(params);
df<-summary(hgOver);
res<-cbind(df[,1],df$Pvalue,df$Count,df$Size,df$Term);
newc<-rep("",dim(res)[1]);
res<-cbind(res,newc);

for(i in 1:dim(res)[1])
{
# for each row, identify the gene names and bind them into the row
gids<- unlist(mget( as.character(geneIdsByCategory(hgOver)[[res[i,1]]] ),envir=org.Hs.egSYMBOL)); 
res[i,6]<-paste(gids,collapse=",");
}
colnames(res)<-c("GOid", "Pvalue", "Count", "Size" ,"Term" ,"Genes");
return(res)
}
################### End of function ####################
########################################################


##########################################################
#   PROCESSING
##########################################################
# take a vector from the frame
universe<-uni[,1]

# Retrieve up and down regulated genes
up<-subset(data,LogFc>0,select=EntrezID)
down<-subset(data,LogFc<0,select=EntrezID)

# If the gene list is in a separate column, 
# you just read this column

#########################################################
## Compute GOenrichment, provide ontology, 
# universe and gene list
# Here ontology biological processes BP, universe in in
# variable universe, gene list is in up
# Other ontologies : Molecular Function MF, Cellular Component CC
# 
#########################################################
rezGplusBP<-GOenrichment("BP",universe,up);

# down
rezGplusBPdown<-GOenrichment("BP",universe,down);

################
# Identify Columns
line<-c("GOid", "Pvalue", "Count", "Size" ,"Term" ,"Genes");
colnames(line)<-c("GOid", "Pvalue", "Count", "Size" ,"Term" ,"Genes");

# Collate all results into one table, separated by line
all<-rbind(rezGplusBP,line,rezGplusBPdown)
write.table(all,file="up_down_GOenrichment.csv",sep="\t",row.names=FALSE)



