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
