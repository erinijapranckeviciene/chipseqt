
###############################################################################################
###  Function scripts that have to be loaded into workspace
###############################################################################################

library("Biostrings")
library("stringr")
library("motifRG")
library("seqLogo")
###############################################################################################
#     HELPER FUNCTIONS FOR STRING MANIPULATION
###############################################################################################
#   1.  Expand consensus string
###############################################################################################
#  expandConsstr <- function (consstr=NULL)
###############################################################################################
#
#  Generates all variants of a consensus sequence fairly quickly. Do not check for mistakes in letters
#  for example Z in the input would stop the script 
#
#  call : r<-expandConsstr(consesnsus_string)
# Parameters: consensus_string= DNA motif consensus string, for example "CAN"
# Output:     string array , for above it would be
#             >r
#              [1] "CAA" "CAC" "CAG" "CAT"
###############################################################################################
expandConsstr <- function (consstr=NULL)
{
    # get IUPAC codes for all letters in a consensus string
    s<-IUPAC_CODE_MAP[unlist(str_split( consstr,pattern=""))][-1]
    # we know that the first letter s[1] is NA
    # initialize sequences by the second letter, remove the first empty element
    sequences<-unlist(str_split(s[1],pattern=""))[-1]
    s<-s[-1];
    for(i in 1:length(s))
    {
        # split the IUPAC_CODES into letters and append them to the sequences
    	letters<-unlist(str_split(s[i],pattern=""))[-1]
    	new_seq<-c();
    	for (j in 1:length(letters)){ new_seq<-c(new_seq, str_pad(sequences,i+1,side="right",letters[j]))}
        sequences<-new_seq;
    }
    
        return(sequences) 
    
}
###############################################################################################
#   2. Transform consensus to regular expressions
###############################################################################################
##  consensusTOregexp<-function(consstr="WGATAAN")
###############################################################################################
#
#   consensusTOregrexp:  converts the consensus string into a regular expression for matching. Removes Ns automtically from the beginnig and end   
#                 
#   Call:  ebox<-consensusTOregexp(consstr="CANNTG")             
#              
#  Parameters: consstr  -  consensus string, for example "CANNTG"
#
#  Output:     regular expression, for above input it is  "CA[A,C,G,T][A,C,G,T]TG"
###############################################################################################

consensusTOregexp<-function(consstr=NULL,trim=TRUE)
{
	input_str<-consstr;
	if(trim){ input_str<-str_trim(str_replace(str_replace(input_str,"^[N]*","  "),"[N]*$","  "),side="both");}
	
	new_s1<-c()
	#split the string into the letters
	s1<-unlist(str_split(input_str,"")[[1]])[1:str_count(input_str,".")+1];
	#print(s1)
	for( i in 1:length(s1))
	{
		letter<-IUPAC_CODE_MAP[s1[i]];
		if( str_count(letter,".")==1){new_s1<-c(new_s1,letter);}
		else
		{
		   new_letter<-str_replace_all(letter,"",",");
		   str_sub(new_letter,start=1,end=1)<-"[";
		   str_sub(new_letter,start=-1,end=-1)<-"]"; 
		   new_s1<-c(new_s1,new_letter);}
	#print(new_s1)	   
	}
	output_str<-str_c(new_s1,collapse="");
	return(output_str)
}
###############################################################################################
#   3. Generate original and reverse patterns and their class as data.frame
#############################################################################################
#          genPatclass<-function(consstr)
#############################################################################################
#genPatclass("WGATAA")
#     pat      class
#[1,] "AGATAA" "0"  
#[2,] "TGATAA" "0"  
#[3,] "AATAGA" "1"  
#[4,] "AATAGT" "1"  
############################################################################################## 
genPatclass<-function(consstr)
{
	pat<-expandConsstr(consstr); class<-rep(0,length(pat)); p1<-cbind(pat,class);
	pat<-expandConsstr(reverse(consstr)); class<-rep(1,length(pat)); p2<-cbind(pat,class);
	return( data.frame(rbind(p1,p2)))
}
###############################################################################################
#   4.  Generate regular expressions for composite motifs at a given distnce 
###############################################################################################
## generateComposites<-function(dist=dist,pat1="CA..TG",pat2="[A,T]GATAA" )
##############################################################################################
#
#   generate_patterns:  given the vector of distances, generates patterns of composite motifs
#                 
#   call:  pat<-generate_patterns(dist=di,pat1=pat1,pat2=pat2)            
#              
#  Parameters: dist   -  vector of distances, for example dist<-seq(from=-20, to=20, by=1);
#              pat1   -  regexp pattern of the first motif 
#              pat2   -  regexp pattern of the secon motif. For negative distances 
#                       the second motif is placed upstream of the first, for positive- downstream  
#  Output:    data.frame of patterns with two columns dist , pattern 
#  >pat
#   dist                                              pattern
#1   -20 [A,T]GATAA....................CA[A,C,G,T][A,C,G,T]TG
#2   -19  [A,T]GATAA...................CA[A,C,G,T][A,C,G,T]TG
# ...
###############################################################################################
generateComposites<-function(dist=seq(from=-20,to=20,by=1),pat1="CA..TG",pat2="[A,T]GATAA" )
{
	# forward patterns
	fp<-c();
	for( i in 1:length(dist))
	{
		if(dist[i] < 0){ p<-paste( c(pat2,paste(rep(".",abs(dist[i])),collapse=""), pat1),collapse="" ); pattern<-cbind(dist[i],p);  } 
		if(dist[i] > 0){ p<-paste( c(pat1,paste(rep(".",abs(dist[i])),collapse=""), pat2),collapse="" ); pattern<-cbind(dist[i],p); }
		if(dist[i] ==0){ patternl<-paste(c(pat2,pat1),collapse=""); patternr<-paste(c(pat1,pat2),collapse="");
			pattern<-rbind( cbind(-0,patternl), cbind( 0, patternr) ) ;}
	fp<-rbind(fp,pattern)	
		 }
	colnames(fp)<-c("dist","pattern");
	fp<-as.data.frame(fp);
	 
	return(fp);
}
###############################################################################################
#   5.  Extracts sequences from peaks, that match original and reverse pattern. Sequences are in the right orientation. 
###############################################################################################
# matching_pattern_in_peaks<-function(peaks,pattern_o,pattern_r)
###############################################################################################
#
#  matching_pattern_in_peaks  
#                  extracts matching sequences from peaks, patterns are
#                  passed as original and reverse, matching is done on original and reverse complement string of peak. 
#                  this function is needed to generate the composite motif PWMs and seqLogos
#  Call:
#        first make the patterns
#
# pr7<-generateComposites(dist=c(-7), pat1=consensusTOregexp("TGNNCA"),pat2=consensusTOregexp("AATAGW"))
# po7<-generateComposites(dist=c(7), pat1=consensusTOregexp("CANNTG"),pat2=consensusTOregexp("WGATAA"))
#
#   call the function
#   pp<-matching_pattern_in_peaks(peaks,paste(po7$pattern),paste(pr7$pattern))
#   
#   make consensus matrix 
#   cm<-consensusMatrix(DNAStringSet(pp))[1:4,]/length(pp)
#   generate logo
#  seqLogo(cm)   
#
###############################################################################################
matching_pattern_in_peaks<-function(peaks,pattern_o,pattern_r)
{  poo<-c(); # peak original pattern original
   por<-c(); # peak original pattern reverse
   pro<-c(); # peak reverse complement , pattern original
   prr<-c(); # peak reverse complement, pattern reverse 
  print(pattern_o);
  print(pattern_r);
   
   matches<-c();	 
	for(i in 1:length(peaks))
	{ 
		peak_o<-toString(peaks[i]);
		peak_r<-toString(reverseComplement(DNAStringSet(peaks[i])));
		# extract pattern as it is
		if (str_detect(peak_o,pattern_o)) 
		{ poo<-c(poo,i); print(paste(c("poo  ",i ),collapse="") );
			d<-str_extract(peak_o,ignore.case(pattern_o)); matches<-c(matches,d); print(d) } 
	    
        # peak original, pattern reversed, we have to reverse resulting pattern
	    if (str_detect(peak_o,pattern_r)) 
	    { por<-c(por,i); print(paste(c("por  ",i ),collapse="") )
	    	d<-str_extract(peak_o,ignore.case(pattern_r)); matches<-c(matches,reverse(d));print(d)}
	    	
	    # peak reverse complement, pattern original, save as it is	 
	    if (str_detect(peak_r,pattern_o))
	     { pro<-c(pro,i); print(paste(c("pro  ",i ),collapse="") )	
	       d<-str_extract(peak_r,ignore.case(pattern_o)); matches<-c(matches,d); print(d) } 
	       
	     # peak reverse complement, pattern reverse, reverse before saving	    
	    if (str_detect(peak_r,pattern_r)) 
	    { prr<-c(prr,i); print(paste(c("prr  ",i ),collapse="") )
	    	d<-str_extract(peak_r,ignore.case(pattern_r)); matches<-c(matches,reverse(d)); print(d)} 
	    
		}
# previous return				
#return(peaks[unique(c( poo, por, pro, prr))]);
return(matches)
}
###############################################################################################
#   6. Information about the matching pattern closest to the summit in a single peak
###############################################################################################
#                     match_pat<-function(fasta,summit,name,pat)
###############################################################################################
# Identifies matching pattern in a fasta string, that is closest to the peak. The eligible 
# patterns are listed in a string array pat
# 
#
# Call : r<-match_pat(peak_fasta_sequence, peak_summit_position, peak_name, pattern)
#        peak_fasta_sequence = sequence underlying peak region  
#        peak_summit_position = peak summit relative to the peak start
#        peak_name = peak name in the form (chr:start-stop)
#        pattern   = data.frame of patterns: pat column is pattern and class column is indicator
#                    whether the pattern is original or reverse. If pattern is reverse, then it's
#                    start  is from the end. This information is used for computing of the distance.  
#                    pattern is output of function gen_pat("pattern"), for example :
#> gen_pat("CAN")
#  pat class
#1 CAA     0
#2 CAC     0
#3 CAG     0
#4 CAT     0
#5 AAC     1
#6 CAC     1
#7 GAC     1
#8 TAC     1      
#
# Example: 
#          rez<-match_pat(peak_fasta_summit$fasta[1],peak_fasta_summit$summit[1],paste(peak_fasta_summit$name[1]),gen_pat("CAN"))
#
#
# Output is a single row of the form :
#  > rez
#> rez
#                      name    patternno    start     stop    pattern     dist     summit 
#"chr1:101628772-101628889"     "4"         "54"      "56"    "CAT"       "-6"      "60" 
# 
#   If match has not been found . then pattern no is assigned -1 and other fields are zero
#   dist is a distance of the pattern start to the summit , negative means that pattern is 
#   upstream of summit
###############################################################################################

match_pat<-function(fasta,summit,name,pattern)
{
# pat should be pattern and category
# if the pattern is reversed, the start position is the pattern's end position
pat<-paste(pattern$pat);
class<-pattern$class;
# locate all occurrences of the patterns in the pat
d1<-str_locate_all(fasta,ignore.case(pat))	
pos<-c(); dist<-c();
#pattern number is i 	
for (i in 1:length(d1))
{   # if there is a match  
	if(length( d1[[i]])>0) 
	{
		d2<-d1[[i]]
		#print(c(d2[1],d2[2]))
	    # get the distance to the summit
		j=1; while(j<=dim(d2)[1])
		{    
			 start<-d2[j,1];
			 stop<-d2[j,2];
			 # if pattern is reversed, then distance is computed from stop
			 if( class[i]==0) { dd<-(start-summit);}
			 else { dd<-(stop-summit); }
			 # record the pattern match
			 dist<-c(dist, abs(dd));
			 s<-toupper(str_sub(fasta,start,stop)); 
			 pos<-rbind(pos,c(name,i,start,stop,s,dd,summit) );
			 j=j+1; } 
	}	
}

# find the minimum
#print(dist)
#print(pos) 
if( length(pos)==0){ r<-c(name,-1,0,0,0,0,summit);}
else{
k<-which.min(dist);
r<-pos[k,];}
names(r)<-c("name","patternno","start","stop","pattern","dist","summit")
return(r)
}
###############################################################################################
#   7. Info about the matching paterns closest to to the summit in all peaks 
#############################################################################################
#   match_pat_all<-function(peaks=peak_fasta_summit, patterns=pat)
############################################################################################## 
#
#  match_pat_all :  finds occurrences of the motif ( original and reversed ) 
#                   in peaks ( original and complement). Consensus sequence of the motif is expanded 
#                   into the explicit set of variants of that motif and set of reverse variants.
#                   The pattern that was closest to the peak summit is recorded. Matches are
#                   reported separate for the original peak ( rez data clumn names are followed by .x ) 
#                   and complement sequence of the peak. ( rez data column names followed by .y )
#  Parameters:
#                    dataframe containing peak name, fasta sequence, patterns to be matched i.e. gen_pat("WGATAA")
#  Output:                  
#                    dataframe with all information about matching pattern. 
#  Example:
#                   call:  rez<-match_pat_all(peaks=peak_fasta_summit, patterns=gen_pat("CANNTG"));   
#> rez[1:3,]
#                      name patternno.x start.x stop.x pattern.x dist.x summit.x patternno.y start.y stop.y pattern.y dist.y summit.y
#1 chr1:101628772-101628889          23      49     54    GTGCAC     -6       60          10      49     54    CACGTG    -11       60
#2   chr1:10699037-10699224          -1       0      0         0      0      136          -1       0      0         0      0      136
#3 chr1:107961734-107961917          -1       0      0         0      0      100          -1       0      0         0      0      100
#
#    The resulting matrix  of match_pat_all() function can be used for ranking motif variants, and cumulative distance plots.  
# 
#    For cumulative distance plot, the minimum distance should be selected min( abs(dist.x, dist.y))
#    For ranking of the motif variants we need separate function 
############################################################################################## 
match_pat_all<-function(peaks=peak_fasta_summit, patterns=pat)
{
rez1<-c(); rez2<-c();
print(patterns)
# for each peak
for (i in 1:dim(peaks)[1])
{   # original sequence
	fasta1<-paste(peaks$fasta[i]);
	#reverse sequence 
	fasta2<-toString(complement(DNAString(paste(fasta1))));
	
	r1<-match_pat(fasta1,peaks$summit[i],paste(peaks$name[i]),patterns); rez1<-rbind(rez1,r1);
	r2<-match_pat(fasta2,peaks$summit[i],paste(peaks$name[i]),patterns); rez2<-rbind(rez2,r2);
 } 	
# merge the two files 
rez<-data.frame(merge(rez1,rez2, by=c("name")))
rownames(rez)<-seq(from=1,to=dim(peaks)[1],by=1);
return(rez)
}
###############################################################################################
#   8.  Retrieves distance to summit of the motifs, closst to summit. Distance is from the 
#       start of the motif till summit position
############################################################################################## 
# motif_dist<-function(motif_rez=rez)
###############################################################################################
#
# motif_dist :    
#                selects motif patterns from file that are closest to the summit 
#                returns data.frame with columns name, pattern, distance 
# 
# Example:    
#            m_dist<-motif_dist(motif_rez=rez);
# Otput:
#>m_dist[1:6,]
#                      name pattern dist
#1 chr1:113601757-113602058  AGATAA  -46
#2   chr1:11518948-11519098  AGATAA   -4
#3   chr1:12221967-12222178  AGATAA  -14
#4   chr1:12664466-12664676  AATAGT  -23
#5 chr1:147266767-147267007  TGATAA   11
#6 chr1:151821168-151821372  AGATAA   20
# ...  
###############################################################################################
motif_dist<-function(motif_rez=rez)
{
mot_dist<-c()	

new_rez<-subset( motif_rez, ((patternno.x!=-1) | (patternno.y!=-1)), select=c("name","pattern.x","dist.x","pattern.y","dist.y") )
# select cases, in which there is only one pattern
x_rez<-subset(new_rez, ( ( pattern.x!=0)&(pattern.y==0)),select=c("name","pattern.x","dist.x")); 
y_rez<-subset(new_rez, ( ( pattern.y!=0)&(pattern.x==0)),select=c("name","pattern.y","dist.y")); 
colnames(x_rez)<-c("name","pattern","dist");
colnames(y_rez)<-c("name","pattern","dist");

n_mot_dist<-data.frame( rbind(x_rez,y_rez) );


# If x_rez or y_rez are empty, it will not have impact
# if both is empty, it will stop execution

both<-subset(new_rez, ( ( pattern.x!=0)&(pattern.y!=0)) ); 
if (dim(both)[1]>0)
{
# in the subset both we have to identify which variant is closer to summit 
 for( i in 1: dim(both)[1])
 { 
 	dx<-as.numeric(paste(both$dist.x[i]));
 	dy<-as.numeric(paste(both$dist.y[i]));
 	patx<-paste(both$pattern.x[i]);
 	paty<-paste(both$pattern.y[i]);
 	name<-paste(both$name[i]);
 	
    if( abs(dx) < abs(dy)){ pattern<-patx; dist<-dx;}
    	else { pattern<-paty; dist<-dy;}
 mot_dist<-rbind(mot_dist,c(name,pattern,dist));
}
colnames(mot_dist)<-c("name","pattern","dist");
n_mot_dist<-rbind( n_mot_dist, data.frame(mot_dist) );
}

colnames(n_mot_dist)<-c("name","pattern","dist");
return(n_mot_dist)
}
###############################################################################################
#   9.  extract fasta strings from DNAStringSet ogject 
############################################################################################## 
#  extract_strings<-function(fasta=fasta) 
###############################################################################################
# fasta is DNAStringSet object, extract data.drame with two columns: name fasta
#
extract_strings<-function(fasta=fasta)
{peak_fasta<-c();
# in order to have only peak name and sequence we have to process DNAStringSet object
for (i in 1:length(fasta))
{peak_fasta<-rbind(peak_fasta,  c( names(fasta)[i], toString(subseq(fasta[[i]]))))}
colnames(peak_fasta)<-c("name","fasta");
return(data.frame(peak_fasta))
}
###############################################################################################
#   10. counts how many given composite motifs are present in peaks 
###############################################################################################
###############################################################################################
#composites2<-function(d=peak_fasta,pat1=pat1,pat2=pat2,di)
###############################################################################################
#
#   composites2: for a set of given peaks and two motifs, composites2() computes the number of occurrences 
#               of a composite motif among peaks. Composite motif is defined as the two motifs, separated  
#               by the distance - given number of bases. The default distances are from -20 to 20. Negative distance means that motif2 is 
#               upstream of the motif1, and positive means the motif2 . 
# 
##  
#  Parameters: d      - data.frame with two columns:  name  fasta
#              pat1   - motif pattern1 as regular expression , generated by function consensusTOregexp 
#              pat2   - motif pattern2 
#              di     - distances
#  Call: ret_rez<-composites2(d=peak_fasta,pat1=pat1,pat2=pat2,di=seq(from=-20, to=20, by=1) )
#  Output:     data.fame of the following form
# ....
##############################################################################################
composites2<-function(d=peak_fasta,pat1=pat1,pat2=pat2,di=seq(from=-20, to=20, by=1) )
{
# Generate patterns original
pat_o<-generateComposites(dist=di,pat1=consensusTOregexp(pat1),pat2=consensusTOregexp(pat2) )
rownames(pat_o)<-seq(from=1,to=dim(pat_o)[1],by=1)
# Count composites for original 
ret_o<-count_composites(data=d,pattern=pat_o)

# Generate patterns reverse
pat_r<-generateComposites(dist=di,pat1=consensusTOregexp(reverse(pat1)),pat2=consensusTOregexp(reverse(pat2)) )
rn<-rownames(pat_r)[dim(pat_r)[1]:1]
rownames(pat_r)<-rn
# Note: count composites takes into account the absolute distance. The pattern name is returned as match
# Count composites for reverse 
ret_r<-count_composites(data=d,pattern=pat_r);

# bind all matches in peaks together
common<-cbind(ret_o[[2]] , ret_r[[2]]);
# initialize counts

cn<-c();

# count composite matches 
for(i in 1:dim(common)[1])
{
    if( any( !is.na(common[i,])))
    { 
    # get the peak matches, m contains name of the composite
    	n<-common[i,]; 
    	
    	m<- as.numeric(paste( n[!is.na(n)] ));
    	
    	#get the distances
    	dt<-as.numeric( paste(pat_o$dist[m]));
    	
    	#get minimum distance
    	k<-which.min(abs(dt) );
    	
    	#remember the name
    	cn<-c(cn, m[k] );
        }
    }
# run length encoding    
t<-rle(cn[order(cn)]);
# create counts
counts<-as.data.frame(cbind(t$values,t$lengths,t$lengths/dim(d)[1]) );
colnames(counts)<-c("name","count","rel_freq")

# get the patterns
pat<-as.data.frame(cbind(rownames(pat_o),pat_o));
colnames(pat)[1]<-c("name");

# merge patterns and counts, keep all rows
# if NA, assign zero count
m<-merge(pat, counts, by=c("name"), all.x=TRUE, all.y=TRUE,sort=FALSE);
m[is.na(m)]=0;
o<-order(as.numeric(paste(m$name))) 
rez<- as.data.frame( cbind( m[o,], ret_o[[1]],ret_r[[1]]));

#assign meaningful names
rownames(rez)<-rez$name;
colnames(rez)<c( "name", "dist","pattern", "count", "rel_freq", "count_oo", "count_ro", "count_or", "count_rr");

#"count_oo"   peak original, composite original, 
#"count_ro",  peak reverse complement, composite original
#"count_or",  peak original, composite reverse 
#"count_rr",  peal reverse complement, composite reverse

return(rez)
} # end of function composites
###############################################################################################
#   11. does the actual counting of the composites
###############################################################################################
###############################################################################################
## count_composites<-function(data=d,pattern=p)
##############################################################################################
#
#   count_composites:  given the data.frame of sequences/name fasta_sequence/ and data.frame of patterns, the function  
#                       counts how many times each individual pattern occurs in peak sequences. If there are several
#                      occurrences, then only the composite with shorter distance is retained               
#   call:  ret<-count_composites(data=d,pattern=pat)           
#              
#  Parameters: data   -  data.frame of peak/other sequences
#              pat    -  data.frame of patterns, output of th function generate_patterns  
#
#  Output:    list with two elements : ret[[1]] is a data.frame in which rows correspond to  the rows of pattern data.frame, 
#                                      columns counts_lr and counts_rl contain counts how many times each pattern occurred 
#                                      in original and reverse complement sequences respectively
#                                      ret[[2]]  is a data frame in which rows correspond to peak sequences, 
#                                      columns matches_lr and matches_rl contains the name(row number of the pattern data.frame)
#                                      of the composite motif that had a match in the original or reverse complement sequence of the peak. 
#  	                                   If there was no match, then NA is assigned. 
###############################################################################################
count_composites<-function(data=d,pattern=po)
{   
	ret<-list();
	counts_lr<-rep(0,dim(pattern)[1]);
	counts_rl<-rep(0,dim(pattern)[1]);
	
	matches_lr<-c();
	matches_rl<-c();	
		
	n<-dim(data)[1]; 
	# find pattern in the data record-string , select the pattern with the smallest distance
	for( i in 1:n) #
	{

	     name_lr<-find_matching_pattern(data_string=paste(data[i,2]),pattern=pattern);
	     
	     name_rl<-find_matching_pattern(data_string=reverseComplement(DNAString(paste(data[i,2]))),pattern=pattern);
	     
	     
	     matches_lr<-c(matches_lr,name_lr);
	     matches_rl<-c(matches_rl,name_rl);
	     
	     if( !is.na(name_lr)){ print(c("lr",i, name_lr)); counts_lr[name_lr]<-counts_lr[name_lr]+1;}
	     if( !is.na(name_rl)){ print(c("rl",i, name_rl)); counts_rl[name_rl]<-counts_rl[name_rl]+1;}		
	        
     }
       
    ret[[1]]<-cbind(counts_lr,counts_rl); 
    ret[[2]]<-cbind(matches_lr,matches_rl);
    return(ret);
}		
###############################################################################################
#   12. does the actual matching of the motif variants to the paek fasta sequence 
###############################################################################################
###############################################################################################
## find_matching_pattern<-function(data_string="NNNNNNN",pattern=pat)
##############################################################################################
#


################################### Find pattern out of the set of the patterns, that is present in the data string
find_matching_pattern<-function(data_string="NNNNNNN",pattern=pat)
{   
	p<-pattern$pattern; dist<-pattern$dist;
    di<-c();ni<-c();
	for(i in 1:length(p))
	{
    	pm<-str_locate_all(toString(data_string),ignore.case(toString(p[i])));
        if(length(pm[[1]])>1) { di<-c(di,abs(as.numeric(paste(dist[i]))) ); ni<-c(ni, i); }
	}
	if(length(di)==0){ return(NA);}
	else {    
		      #return(  ni[which.min(di)] ); 
		      element<- ni[which.min(di)];
		      # return the name of the pattern. 
		      return( as.numeric(rownames(pattern)[element]) )
		      
		      
		      }
}
###############################################################################################
#   13. plot preferred distances
###############################################################################################
#####################   Preferred distances Plot function ################################
## plot_pd<-function(dist,peak_pd,bg_pd,tit=tit,maxy=maxy,name=name)
#  dist - distances on x axis, dist= dist<-bg[,1];
#  peak_pd - relative frequency of the composite in peaks
#  bg_pd   - relative frequency of the composite in background
#  tit -     title on the plot
#  name-     name on the plot and a file name to save pdf
######################################################################
plot_pd<-function(dist,peak_pd,bg_pd,tit=tit,maxy=maxy,name=name)
{ 

# Plot for peak and background
# Plots are in png files 
pdf(eval( paste(c(name,".pdf"),collapse="") ),title="TAL1 ECFC",paper="special",width=12,height=7)
par(mfrow=c(1,2))
plot( dist, peak_pd,type="h",xlim=c(-20,20), ylim=c(0,maxy),xlab="Distance (bp)", ylab="Relative Frequency", main=tit); grid();
plot( dist, bg_pd,type="h",xlim=c(-20,20), ylim=c(0,maxy),xlab="Distance (bp)", ylab="Relative Frequency", main="background"); grid();
dev.off()
}
###############################################################################################
#   14.   plot the distances from the summit of the motifs that are closest from the summit 
###############################################################################################
#  cum_dist_plot<-function(di=dist,tot=tot,m_name=m_name)
###############################################################################################

# cum_dist_plot  : plots the distances as cumulative distribution function into pdf file
#                  
#   Parameters:    di = vector of distances
#                  tot= total number of peaks
#                  m_name= motif name for plot and file name
################################################################################################
cum_dist_plot<-function(di=dist,tot=tot,m_name=m_name)
{
   # create run_length encoded variable  
   # NOTE: data.frame variables usualy are factors. We have to transform them into numeric
     di<-abs(as.numeric(paste(di)));   	
     o<-order(di);
     r<-rle(di[o]);
     print(r)
     values<-as.numeric(paste(r$values))
     lengths<-as.numeric(paste(r$lengths))
     cdf<-c(cumsum(lengths)/tot);
pdf(eval( paste(c(m_name,".pdf"),collapse="") ),title="   motif distance to the peak summit")
plot(values,cdf,type='b',col='blue',ylim=c(0,1),xlim=c(0,300),xlab="Distance bp", ylab=" Total fraction of peaks",lwd=5,main=m_name )
dev.off()
}    

###############################################################################################
#   15.   select only those peaks, that have a pattern, peaks are in DNAStringSet 
###############################################################################################
# matching_peaks<-function(peaks,pattern)
###############################################################################################
#
#  matching_peaks  tests for the presence of the motif pattern in
#                  original and reverse complement string of peak. It
#                  returns only peaks that have matches. This function is called by motifPWMs
#  Parameters:     peaks -  DNAStringSet of peaks
#                  pattern - regular expression pattern of the consensus string  

#  Call:  r<-matching_peaks(peaks,consensusTOregexp("WGATAA"))
#  Output: r = matching peaks in DNAStringSet format  
###############################################################################################
matching_peaks<-function(peaks,pattern)
{  nf<-c(); 
	for(i in 1:length(peaks))
	{ if( (str_detect(toString(peaks[i]),pattern)) || (str_detect( toString(reverseComplement(peaks[i])),pattern) )  ){ nf<-c(nf,i); print(i)} 
	  
		}
return(peaks[unique(c(nf))]);}
###############################################################################################


###############################################################################################
#   16.   Generate PWM model for de novo/ or known motifs by Zihzen's package
###############################################################################################
# motifPWMs<-function(mseq,peaks)
###############################################################################################
#
#  motifPWMs  generates position weight matrix for consensus sequences , using the original peaks 
#             it returns a list of two elements > r[[1]] contains text summarizing all motifs,
#             r[[2]] contains list of PWMs in order to use them in generating  Sequence Logos 
#             This funtion refers to a  function in Zizhen's motifRG package , which requires that the consensus 
#             sequence contain all letters (A,C.G,T). If, upon expansion of the consensus
#             sequence not all four letter are present, in order to avoid interruption of batch processing
#             the array of expanded consensus sequences is augmented by N =[ A,C.G,T ] . For example WGATAA
#             is expanded into "AGATAA" "TGATAA". The C is missing. We append the N to the end of the sequences   
#             and correct back to original length after computations. PWMs of the motifs are identified using the peaks,
#              that have matches of the motif.
#             
#  Call: r<-motifPWMs( mseq, peaks )
#
#  Parameters:  mseq   = array of consensus motif sequences         
#               peaks  = peaks as DNAStringset
#   
#  Output:      r list of two elements
#  Example: 
#              r<-motifPWMs("WGATAA",peaks);
#  First element of list contains plain text lines, that can be submitted to Jaspar 
#> print(r[[1]])
#[1] "> WGATAA"                       
#[2] "A [ 61 0.2 99.4 0.4 99.4 99.4 ]"
#[3] "C [ 0.2 0.2 0.2 0.2 0.2 0.2 ]"  
#[4] "G [ 0.4 99.4 0.2 0.2 0.2 0.2 ]" 
#[5] "T [ 38.4 0.2 0.2 99.2 0.2 0.2 ]"

##
##  Second element of the list contains PWM matrix, that can be used to generte LOGO
#> print(r[[2]][[1]])
#  [,1] [,2] [,3] [,4] [,5] [,6]
#A 61.0  0.2 99.4  0.4 99.4 99.4
#C  0.2  0.2  0.2  0.2  0.2  0.2
#G  0.4 99.4  0.2  0.2  0.2  0.2
#T 38.4  0.2  0.2 99.2  0.2  0.2

#
#  command to enerate the SeqLogo
# 
#> seqLogo(r[[2]][[1]]/100)
#
########################################################################
motifPWMs<-function(mseq,peaks)
{  rez_all<-c() 
   rez_pwm<-list() 

# since some motifs are underrepresented in peaks,  then  
# identification of the motif PWM by using the peaks that do not conain 
# that motif lead to incorrect representation of position weights. 
# for this reason we use only peaks that contain a match to the consensus sequence      
#  

   for (i in 1:length(mseq))
	{
        new_peaks<-matching_peaks(peaks,consensusTOregexp(mseq[i]));
		p<-expandConsstr(mseq[i])	
		correct<-FALSE;
		# check whether consensus string has all four letters
		if (0 %in% str_count(paste(p,collapse=""),c("A","C","G","T"))) {
			new_p<-str_pad(p,str_length(p[1])+1,side="right",c("A"));
			new_p<-c(new_p, str_pad(p,str_length(p[1])+1,side="right",c("C")));
			new_p<-c(new_p, str_pad(p,str_length(p[1])+1,side="right",c("G")));
            new_p<-c(new_p, str_pad(p,str_length(p[1])+1,side="right",c("T")));
            p<-new_p; correct<-TRUE;
			};			
		pm<-refinePWMMotif(p,new_peaks,max.iter=500)
		# get the PWM rounded to two digits
		rez<-round(pm$model$prob*100,digits=2);
		if(correct){ no<-dim(rez)[2]; rez<-rez[,-no];}
		# save PWM
		rez_pwm[[i]]<-rez;
        # generate the text
		name<-paste(c(">",mseq[i]),collapse=" ")
		rez_all<-c(rez_all,name)
        for (k in 1:dim(rez)[1])
		 { 
		 	o<-paste(rez[k,]); oo<-paste(c(rownames(rez)[k],"[",o,"]"),collapse=" ")
		 	rez_all<-c(rez_all,oo)
          }
    }
ret_list<-list();
ret_list[[1]]<-rez_all;
ret_list[[2]]<-rez_pwm;    
return(ret_list)
}
#############################################################################################
### NEW FUNCTIONS
##############################################################################################
#   10. counts how many given composite motifs are present in peaks 
###############################################################################################
###############################################################################################
#composites2<-function(d=peak_fasta,pat1=pat1,pat2=pat2,di)
###############################################################################################
#
#   composites2: for a set of given peaks and two motifs, composites2() computes the number of occurrences 
#               of a composite motif among peaks. Composite motif is defined as the two motifs, separated  
#               by the distance - given number of bases. The default distances are from -20 to 20. Negative distance means that motif2 is 
#               upstream of the motif1, and positive means the motif2 . 
# 
##  
#  Parameters: d      - data.frame with two columns:  name  fasta
#              pat1   - motif pattern1 as regular expression , generated by function consensusTOregexp 
#              pat2   - motif pattern2 
#              di     - distances
#  Call: ret_rez<-composites2(d=peak_fasta,pat1=pat1,pat2=pat2,di=seq(from=-20, to=20, by=1) )
#  Output:     data.fame of the following form
# ....
##############################################################################################
composites3<-function(d=peak_fasta,pat1=pat1,pat2=pat2,di=seq(from=-20, to=20, by=1) )
{
# Generate patterns original
pat_o<-generateComposites(dist=di,pat1=consensusTOregexp(pat1),pat2=consensusTOregexp(pat2) )
rownames(pat_o)<-seq(from=1,to=dim(pat_o)[1],by=1)
# Count composites for original 
ret_o<-count_composites(data=d,pattern=pat_o)

# Generate patterns reverse
pat_r<-generateComposites(dist=di,pat1=consensusTOregexp(reverse(pat1)),pat2=consensusTOregexp(reverse(pat2)) )
rn<-rownames(pat_r)[dim(pat_r)[1]:1]
rownames(pat_r)<-rn
# Note: count composites takes into account the absolute distance. The pattern name is returned as match
# Count composites for reverse 
ret_r<-count_composites(data=d,pattern=pat_r);

# bind all matches in peaks together
common<-cbind(ret_o[[2]] , ret_r[[2]]);
# initialize counts

cn<-c();

# count composite matches 
for(i in 1:dim(common)[1])
{
    if( any( !is.na(common[i,])))
    { 
    # get the peak matches, m contains name of the composite
    	n<-common[i,]; 
    	
    	m<- as.numeric(paste( n[!is.na(n)] ));
    	
    	#get the distances
    	dt<-as.numeric( paste(pat_o$dist[m]));
    	
    	#get minimum distance
    	k<-which.min(abs(dt) );
    	
    	#remember the name
    	cn<-c(cn, m[k] );
        }
    }
# run length encoding    
t<-rle(cn[order(cn)]);
# create counts
counts<-as.data.frame(cbind(t$values,t$lengths,t$lengths/dim(d)[1]) );
colnames(counts)<-c("name","count","rel_freq")

# get the patterns
pat<-as.data.frame(cbind(rownames(pat_o),pat_o));
colnames(pat)[1]<-c("name");

# merge patterns and counts, keep all rows
# if NA, assign zero count
m<-merge(pat, counts, by=c("name"), all.x=TRUE, all.y=TRUE,sort=FALSE);
m[is.na(m)]=0;
o<-order(as.numeric(paste(m$name))) 
rez<- as.data.frame( cbind( m[o,], ret_o[[1]],ret_r[[1]]));

#assign meaningful names
rownames(rez)<-rez$name;
colnames(rez)<c( "name", "dist","pattern", "count", "rel_freq", "count_oo", "count_ro", "count_or", "count_rr");

#"count_oo"   peak original, composite original, 
#"count_ro",  peak reverse complement, composite original
#"count_or",  peak original, composite reverse 
#"count_rr",  peal reverse complement, composite reverse

return(rez)
} # end of function composites
######

#############
#   11. does the actual counting of the composites
###############################################################################################
###############################################################################################
## count_composites3<-function(data=d,pattern=p)
##############################################################################################
#
#   count_composites:  given the data.frame of sequences/name fasta_sequence/ and data.frame of patterns, the function  
#                       counts how many times each individual pattern occurs in peak sequences. If there are several
#                      occurrences, then only the composite with shorter distance is retained               
#   call:  ret<-count_composites3(data=d,pattern=pat)           
#              
#  Parameters: data   -  data.frame of peak/other sequences
#              pat    -  data.frame of patterns, output of th function generate_patterns  
#
#  we need to count all occurrences of the patterns in all variants of the string
#  we exhaustively register everything
#
###############################################################################################
count_composites3<-function(data=d,pat=po)
{   
	# will create data frame: pattern, distance, position, orientation, peak name 
	#
	pmatches<-c();
	#
	n<-dim(data)[1]; 
	# find pattern in the data_string, do not select the pattern with the smallest distance, register all
	for( i in 1:n) #all peaks 
	{    data_string=paste(data[i,2]);
         data_name=paste(data[i,1]);
         
##################### find_matching_pattern

	      p<-pat$pattern; dist<-pat$dist;
          # register each pattern if the match is found
          # pattern goes from left to the right and from the right to the left on original and complement strands
          original_lr<-toString(data_string);
          original_rl<-toString(reverse(data_string));
          complement_lr<-toString(complement(DNAString(data_string)));
           complement_rl<-toString(reverse(complement_lr)); 
         
	      
	      for(j in 1:length(p))
	      {
    	   # match on original_lr 
    	   #print(j)
    	   pp=toString(p[j]);
    	   pm1<-str_locate_all(original_lr,ignore.case(pp));
    	   if(length(pm1[[1]])>0)
    	   {
     	   for(k in 1:length(pm1))
    	      { st=pm1[[k]][1]; en=pm1[[k]][2]; row<-c(i,data_name,"o_lr",st,en,str_length(original_lr),paste(dist[j]),pp);       
    	   	           print(paste(row,collapse=" "));
    	   	           pmatches<-rbind(pmatches,row);}
    	   }
    	   # match on original_rl
    	
    	   pm2<-str_locate_all(original_rl,ignore.case(toString(p[j])));
          if(length(pm2[[1]])>0)
    	   {
    	   for(k in 1:length(pm2))
    	      { st=pm2[[k]][1]; en=pm2[[k]][2]; row<-c(i,data_name,"o_rl",str_length(original_rl)-en,str_length(original_rl)-st,str_length(original_rl),paste(dist[j]),pp);       
    	   	           print(paste(row,collapse=" "));
    	   	           pmatches<-rbind(pmatches,row);}
    	   }
 
    	   # match on complement_lr
    	   pm3<-str_locate_all(complement_lr,ignore.case(toString(p[j])));
          if(length(pm3[[1]])>0)
    	   {
    	    
    	    for(k in 1:length(pm3))
    	      { st=pm3[[k]][1]; en=pm3[[k]][2]; row<-c(i,data_name,"c_lr",st,en,str_length(complement_lr),paste(dist[j]),pp);       
    	   	           print(paste(row,collapse=" "));
    	   	           pmatches<-rbind(pmatches,row);}
    	    
    	    	   }

           # match on complement_rl
    	   pm4<-str_locate_all(complement_rl,ignore.case(toString(p[j])));
          if(length(pm4[[1]])>0)
    	   {
             for(k in 1:length(pm4))
    	       { st=pm4[[k]][1]; en=pm4[[k]][2]; row<-c(i,data_name,"c_rl",str_length(complement_rl)-en,str_length(complement_rl)-st,str_length(complement_rl),paste(dist[j]),pp);       
    	   	           print(paste(row,collapse=" "));
    	   	           pmatches<-rbind(pmatches,row);}

    	   }

          }
}
pm<-as.data.frame(pmatches);
colnames(pm)<-c("peak_no","peak_name","orient","mstart","mstop","peak_len","dist","pattern");	
rownames(pm)<-seq(from=1, to=dim(pm)[1]);

########################################
ret<-list();	     
    
n<-as.numeric(paste(pm$dist))
o=order(n)
nrle<-rle(n[o])

ret[[1]]<-pm;
ret[[2]]<-nrle;
#plot(nrle$values,nrle$lengths,type="h") 
    
    
    return(ret);
}		
######

#############
#   12. does the actual matching of the motif variants to the paek fasta sequence 
###############################################################################################
###############################################################################################
## find_matching_pattern<-function(data_string="NNNNNNN",pattern=pat)
##############################################################################################
#


################################### Find pattern out of the set of the patterns, that is present in the data string
find_matching_pattern3<-function(data_string="NNNNNNN",pattern=pat)
{   
	p<-pattern$pattern; dist<-pattern$dist;
    di<-c();ni<-c();
	for(i in 1:length(p))
	{
    	pm<-str_locate_all(toString(data_string),ignore.case(toString(p[i])));
        if(length(pm[[1]])>1) { di<-c(di,abs(as.numeric(paste(dist[i]))) ); ni<-c(ni, i); }
	}
	if(length(di)==0){ return(NA);}
	else {    
		      #return(  ni[which.min(di)] ); 
		      element<- ni[which.min(di)];
		      # return the name of the pattern. 
		      return( as.numeric(rownames(pattern)[element]) )
		      
		      
		      }
}
###
###############################################################################################
pattern_in_peaks<-function(peaks,pattern_o,pattern_r)
{  poo<-c(); # peak original pattern original
   por<-c(); # peak original pattern reverse
   pro<-c(); # peak reverse complement , pattern original
   prr<-c(); # peak reverse complement, pattern reverse 
  print(pattern_o);
  print(pattern_r);
   
   matches<-c();	 
	for(i in 1:length(peaks))
	{ 
		peak_o<-toString(peaks[i]);
		peak_r<-toString(reverse(DNAStringSet(peaks[i])));
		# extract pattern as it is
		if (str_detect(peak_o,pattern_o)) 
		{ poo<-c(poo,i); print(paste(c("poo  ",i ),collapse="") );
			d<-str_extract(peak_o,ignore.case(pattern_o)); matches<-c(matches,d); print(d) } 
	    
        # peak original, pattern reversed, we have to reverse resulting pattern
	    if (str_detect(peak_o,pattern_r)) 
	    { por<-c(por,i); print(paste(c("por  ",i ),collapse="") )
	    	d<-str_extract(peak_o,ignore.case(pattern_r)); matches<-c(matches,reverse(d));print(reverse(d))}
	    	
	    # peak reverse complement, pattern original, save as it is	 
	    if (str_detect(peak_r,pattern_o))
	     { pro<-c(pro,i); print(paste(c("pro  ",i ),collapse="") )	
	       d<-str_extract(peak_r,ignore.case(pattern_o)); matches<-c(matches,d); print(d) } 
	       
	     # peak reverse complement, pattern reverse, reverse before saving	    
	    if (str_detect(peak_r,pattern_r)) 
	    { prr<-c(prr,i); print(paste(c("prr  ",i ),collapse="") )
	    	d<-str_extract(peak_r,ignore.case(pattern_r)); matches<-c(matches,reverse(d)); print(reverse(d))} 
	    
		}
# previous return				
#return(peaks[unique(c( poo, por, pro, prr))]);
return(matches)
}

####################################################
# LOCATE motif patterns in the string
####################################################
locate_patterns<-function(st,forw,rev)
####################################################
#  st   - string DNA letters
#  forw - short sequence 
#  rev  - reversed forw sequence
# the sequences are regular expressions
####################################################
{

#########################
start_end<-function(p)
#########################
{
ss<-c()	
u=unlist(p); n=length(u); step=(n/2); st=seq(1:step); en=st+step; start=u[st]; end=u[en];
for(i in 1:length(start)){ ss<-c(ss,paste(start[i],end[i],collapse=" ") ) }
r=paste(ss,collapse=",")
return(r)
}

####################################
get_patterns<-function(st,forw,cls)
####################################
{
    patterns  <-c()
    positions <-c()
    class     <-c()    
     r<-list()
     
	if(str_detect(st,ignore.case(forw)))
	{  #get the locations 
	   pl<-str_locate_all(st,ignore.case(forw))
	   #get the patterns
	   pp<-unlist(str_extract_all(st,ignore.case(forw))) 
       s<-start_end(pl)
       c=rep(cls,length(pp))
    
    r[[1]]<-pp
    r[[2]]<-s
    r[[3]]<-c
    
    return(r)
    
	}
}

############# THE COMPLEMENT OF THE PEAK############
    stc<-toString(complement(DNAString(st))) 

############# VARIABLES ############
patterns  <-c()
positions <-c()
class     <-c()    
    
	# CLASS ff
	# DNA string forward
	# The pattern forward
    rff<-get_patterns(st,forw,"ff")
    
    
	# CLASS fr
	# DNA string forward
	# The pattern reverse

    rfr<-get_patterns(st,rev,"fr")
    
	# CLASS cf
	# DNA string complement 
	# The pattern forward
    rcf<-get_patterns(stc,forw,"cf")


	# CLASS cr
	# DNA string complement
	# The pattern reverse
    rcr<-get_patterns(stc,rev,"cf")


if(length(rff)>0)
{ patterns<-c(patterns,rff[[1]])
  positions<-c(positions,rff[[2]])
  class<-c(class,rff[[3]])
	}

if(length(rfr)>0)
{ patterns<-c(patterns,rfr[[1]])
  positions<-c(positions,rfr[[2]])
  class<-c(class,rfr[[3]])
	}	
	
	
if(length(rcf)>0)
{ patterns<-c(patterns,rcf[[1]])
  positions<-c(positions,rcf[[2]])
  class<-c(class,rcf[[3]])
	}

if(length(rcr)>0)
{ patterns<-c(patterns,rcr[[1]])
  positions<-c(positions,rcr[[2]])
  class<-c(class,rcr[[3]])
	}	
	
r<-list()
r[[1]]<-patterns
r[[2]]<-positions
r[[3]]<-class	

return(r)	
	
	}
##################################################
patterns_in_peaks<-function(peaks,forwstr,revstr)
##################################################
{
	n=dim(peaks)[1]
    d<-matrix(nrow=n,ncol=4)
	
    for(i in 1:n)
     {   peakstr<-toString(peaks[i,2])
     	peakname<-peaks[i,1]
     	
     	p<-locate_patterns(peakstr,forwstr,revstr)
     	
     	if(length(p)>0)
     	{
     	patterns<-paste(p[[1]],collapse=",")
     	positions<-paste(p[[2]],collapse=",")
     	class<-paste(p[[3]],collapse=",")
     	d[i,1]=paste(peakname)
     	d[i,2]=patterns
     	d[i,3]=positions
     	d[i,4]=class
         }
         else{
         d[i,1]=paste(peakname)
     	d[i,2]="NA"
     	d[i,3]="NA"
     	d[i,4]="NA"	
         	
         }
         print(i)
	}
colnames(d)<-c("name","patterns","positions","class")
return(d)	
}

