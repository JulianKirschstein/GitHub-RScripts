#plots for foldchanges
#Hallo
##

foldchangebetweentworuns=function(untreated1,untreated2,treated1,treated2){     # new function that requires 4 datasets, 2 untreated and 2 treated (duplicates)
                  
  genecountuntreated1=c(untreated1$match)               #
  genecountuntreated2=c(untreated2$match)               # new vectors which contain the number of matches
  genecounttreated1=c(treated1$match)                   # for the 4 different datasets
  genecounttreated2=c(treated2$match)                   #
                          
  sumuntreated1=sum(untreated1$match)                   #
  sumuntreated2=sum(untreated2$match)                   # compute the sum of all matches in one run which is
  sumtreated1=sum(treated1$match)                       # needed for normalisation between different conditions
  sumtreated2=sum(treated2$match)                       #
                          
  sumuntreated=sumuntreated1+sumuntreated2              # combine the sums of runs with the same conditions
  sumtreated=sumtreated1+sumtreated2                    #
  
  genecountuntreatedmean=c(NA,NA)
  genecounttreatedmean=c(NA,NA)
  
  for(i in 1:length(genecountuntreated1)){                              #
    bothgenecounts1=c(genecountuntreated1[i],genecountuntreated2[i])    # compute the mean genecount (match) for 
    genecountuntreatedmean[i]=mean(bothgenecounts1)                     # runs with condition "untreated"
                            
    if(genecountuntreatedmean[i]<0.5){                                  # set value of mean genecount !=0
      genecountuntreatedmean[i]=0.0001                                  # (division by 0 causes problems)
    }                                                                   
  }                                                                   
                          
  for(i in 1:length(genecounttreated1)){                                #                           
    bothgenecounts2=c(genecounttreated1[i],genecounttreated2[i])        # compute the mean genecount (match) for
    genecounttreatedmean[i]=mean(bothgenecounts2)                       # runs with condition "treated"
                                                                
    if(genecounttreatedmean[i]<0.5){                                    # set value of mean genecount !=0
      genecounttreatedmean[i]=0.0001                                    # (division by 0 causes problems)    
    }                                     
  }                                                                     
                          
  genecountuntreatedmeannorm=genecountuntreatedmean/sumuntreated        # normalisation of mean genecounts             
  genecounttreatedmeannorm=genecounttreatedmean/sumtreated              # (division by total genecount)
  
  foldchange=genecounttreatedmeannorm/genecountuntreatedmeannorm        # 
  log10foldchange=log10(foldchange)                                     # compute foldchange,log10(foldchange) and log2(foldchange)
  log2foldchange=log2(foldchange)                                       #
                          
  cols=c("red","green")                                                                                                             #
  pos=log10foldchange>=0                                                                                                            # plot log10(foldchange vs. genes)
  barplot(log10foldchange,ylim=c(-5,2),border=cols[pos+1],xlab="genes",ylab="fold change (log10)",main="difference (fold change)")  #
                          
  cols=c("red","green")                                                                                                             #
  pos=log2foldchange>=0                                                                                                             # plot log2(foldchange vs. genes)
  barplot(log2foldchange,ylim=c(-15,5),border=cols[pos+1],xlab="genes",ylab="fold change (log2)",main="difference (fold change)")   #
}

mean.foldchange=function(untreated1,treated1){              # new function that requires 2 datasets (untreated and treated)
  
  genecount=c(untreated1$match)                             # new vectors which contain data
  designspergene=c(untreated1$X..designs.total)             # to work with
  meanpergene=c(NA,NA)                                                        
  
  for(i in 1:length(genecount)){                            # compute mean count per design for  
    meanpergene[i]=genecount[i]/designspergene[i]           # every gene under "untreated" condition
    }
                
  genecountuntreated1=c(untreated1$match)                   # vectors with number of read for all genes
  genecounttreated1=c(treated1$match)                       # ("untreated" and "treated" condition)
  
  sumuntreated1=sum(untreated1$match)                       # total number of reads per condition
  sumtreated1=sum(treated1$match)                           # (needed for normalisation)
                
  for(i in 1:length(genecountuntreated1)){                  #
    if(genecountuntreated1[i]<0.5){                         # set value of genecount !=0
      genecountuntreated1[i]=0.0001                         # (division by 0 causes problems)
    }
  }
                
  for(i in 1:length(genecounttreated1)){                    #
    if(genecounttreated1[i]<0.5){                           # set value of genecount !=0
      genecounttreated1[i]=0.0001                           # (division by 0 causes problems)
    }
  }
  
  genecountuntreated1norm=genecountuntreated1/sumuntreated1   # normalistaion of genecount
  genecounttreated1norm=genecounttreated1/sumtreated1         # (division by total genecount)
  
  foldchange=genecounttreated1norm/genecountuntreated1norm    # compute foldchange and log2(foldchange)
  log2foldchange=log2(foldchange)                             #
                
  plot(meanpergene,log2foldchange,pch=20,col="red",ylab="log2 fold change",xlab="mean (per design per gene)",main="differential expression")   # plot log2(foldchange vs. mean per design per gene)
}

foldchange.coverage=function(untreated1,treated1){      # new function that requires 2 datasets (untreated and treated)
                    
  readcount=c(untreated1$match)                         #
  genecountuntreated1=c(untreated1$match)               # vectors of data to work with
  genecounttreated1=c(treated1$match)                   # 
  
  sumuntreated1=sum(untreated1$match)                   # total number of reads per condition
  sumtreated1=sum(treated1$match)                       # (needed for normalisation)
                    
  for(i in 1:length(genecountuntreated1)){              #
    if(genecountuntreated1[i]<0.5){                     # set value of genecount !=0
      genecountuntreated1[i]=0.0001                     # (division by 0 causes problems)
    }
  }
                    
  for(i in 1:length(genecounttreated1)){                #
    if(genecounttreated1[i]<0.5){                       # set value of genecount !=0
      genecounttreated1[i]=0.0001                       # (division by 0 causes problems)
    }
  }
                    
  genecountuntreated1norm=genecountuntreated1/sumuntreated1     # normalistaion of genecount
  genecounttreated1norm=genecounttreated1/sumtreated1           # (division by 0 causes problems)
  
  foldchange=genecounttreated1norm/genecountuntreated1norm      # compute foldchange and log2(foldchange)
  log2foldchange=log2(foldchange)                               #

  readcount.log2foldchange=c(readcount,log2foldchange)              # matrix (nx2) with
  dim(readcount.log2foldchange)=c(length(readcount),2)              # readcounts and log2(foldchange)
  
  readcount.log2foldchange=as.data.frame(readcount.log2foldchange)  # transform matrix 
  names(readcount.log2foldchange)=c("reads","foldchange")           # into data.frame
  
  order.reads=order(readcount.log2foldchange$reads)                       # sort by readcounts
  
  sortedreadcount.log2foldchange=readcount.log2foldchange[order.reads,]   # new data.frame with sorted redcounts
  
  sortedfoldchange=(sortedreadcount.log2foldchange$foldchange)            # vector with sorted foldchanges
  sortedreadcount=(sortedreadcount.log2foldchange$reads)                  # vector with sorted readcounts
  
  df.frequenciesofsamereadcounts=as.data.frame(table(sortedreadcount))    # data.frame with readcounts and their frequency
  names(df.frequenciesofsamereadcounts)=c("readcounts","frequency")       #
  #print(readcount.log2foldchange)
  #print(df.frequenciesofsamereadcounts)
  numberofdifferentreadcounts=length(df.frequenciesofsamereadcounts$frequency)  # number of DIFFERENT frequencies
            
  allfrequencies=c(df.frequenciesofsamereadcounts$frequency)                    # vector with all different frequncies
  meanfoldchangeforonenumberofreadcount=c(NA,NA)
  
  i=1                                                                                       #
  for(j in 1:numberofdifferentreadcounts){                                                  #
    allfoldchangesforonenumberofreadcount=c(sortedfoldchange[i:(i+allfrequencies[j]-1)])    # compute the mean foldchange for
    meanfoldchangeforonenumberofreadcount[j]=mean(allfoldchangesforonenumberofreadcount)    # all different readcounts
    i=(i+allfrequencies[j])                                                                 #
  }
                    
  #allreadcounts=c(sortedreadcount.log2foldchange$reads)                                            #
  #uniquereadcounts=unique(allreadcounts)                                                           # data.frame with readcounts
  #foldchanges.uniquereadcounts=data.frame(uniquereadcounts,meanfoldchangeforonenumberofreadcount)  # and mean foldchange
  #names(foldchanges.uniquereadcounts)=c("readcount","mean-foldchange")                             #
                    
  cols=c("red","green")                                                                                                                                       #           
  pos=meanfoldchangeforonenumberofreadcount>=0                                                                                                                #
  barplot(meanfoldchangeforonenumberofreadcount,border=cols[pos+1],xlab="number of read",ylab="mean log2 fold change",main="mean fold change vs. coverage")   # plot mean log2(folchanges) vs. number of reads
}

sorted.folchange.of.genes=function(untreated1,treated1){
  
  readcount=c(untreated1$match)                         #
  genecountuntreated1=c(untreated1$match)               # vectors of data to work with
  genecounttreated1=c(treated1$match)                   # 
  gID=c(untreated1$GeneID)
  
  sumuntreated1=sum(untreated1$match)                   # total number of reads per condition
  sumtreated1=sum(treated1$match)                       # (needed for normalisation)
  
  for(i in 1:length(genecountuntreated1)){              #
    if(genecountuntreated1[i]<0.5){                     # set value of genecount !=0
      genecountuntreated1[i]=0.0001                     # (division by 0 causes problems)
    }
  }
  
  for(i in 1:length(genecounttreated1)){                #
    if(genecounttreated1[i]<0.5){                       # set value of genecount !=0
      genecounttreated1[i]=0.0001                       # (division by 0 causes problems)
    }
  }
  
  genecountuntreated1norm=genecountuntreated1/sumuntreated1         # normalistaion of genecount
  genecounttreated1norm=genecounttreated1/sumtreated1               # (division by 0 causes problems)
  
  foldchange=genecounttreated1norm/genecountuntreated1norm          # compute foldchange and log2(foldchange)
  log2foldchange=log2(foldchange)                                   #
  
  geneID.log2foldchange=c(gID,log2foldchange)                       # matrix (nx2) with
  dim(geneID.log2foldchange)=c(length(readcount),2)                 # readcounts and log2(foldchange)
  
  geneID.log2foldchange=as.data.frame(geneID.log2foldchange)        # transform matrix 
  names(geneID.log2foldchange)=c("geneID","foldchange")             # into data.frame
  
  order.foldchanges=order(geneID.log2foldchange$foldchange)         # sort by foldchanges
  
  geneID.sortedlog2foldchange=geneID.log2foldchange[order.foldchanges,]   # new data.frame with sorted foldchanges
  
  orderedfoldchanges=c(geneID.sortedlog2foldchange$foldchange)
  
  plot(1:length(gID),orderedfoldchanges,pch=20,xlab="genes (ordered by foldchange)",ylab="log2 foldchange",main="foldchange (sorted)")
}

foldchangebetweentworuns(untreated1,untreated2,treated1,treated2)  # 
mean.foldchange(untreated1,treated1)                               # execute all the functions 
foldchange.coverage(untreated1,treated1)                           # 
sorted.folchange.of.genes(untreated1,treated1)                     #