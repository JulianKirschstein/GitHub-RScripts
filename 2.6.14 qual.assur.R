library(plotrix)

#load "2132_CRISPR_KO_Plasmid-designs" as "PlasmidDesigns"
#load "2132_CRISPR_KO_Plasmid-genes" as "PlasmidGenes"


#quality assurance: plot number of reads vs. number of targets

read.distr.targets=function(PlasmidDesigns){        # new function that requires 1 dataset  (PlasmidDesigns)
  
  targetcount=c(PlasmidDesigns$match)               # vector with all targetcounts
  #meandesigns=mean(targetcount)                           
  #cat("The mean of reads per target is",meandesigns,"\n")
  
  hist(targetcount,breaks=100,col="blue",xlab="number of reads (per target)",ylab="number of targets",main="read distribution")                   # plot number of targets vs. number of reads (read distribution)
  hist(targetcount,breaks=200,xlim=c(0,3000),col="blue",,xlab="number of reads (per target)",ylab="number of targets",main="read distribution")   # with 2 different limits for the x-Axis
}

#quality assurance: plot number of reads vs. number of genes

read.distr.genes=function(PlasmidGenes){        # new function that requires 1 dataset (PlasmidGenes)
  
  genecount=c(PlasmidGenes$match)               # vector with all genecounts
  #meangenes=mean(genecount)
  #cat("The mean of reads per gene is",meangenes,"\n")
  
  hist(genecount,breaks=100,col="blue",xlab="number of reads (per gene)",ylab="number of genes",main="read distribution")                     # plot number of genes vs. number of reads (read distribution)
  hist(genecount,breaks=1000,xlim=c(0,15000),col="blue",,xlab="number of reads (per gene)",ylab="number of genes",main="read distribution")   # with 2 different limits for the x-Axis
}

#quality assurance: plot number of targets vs. number of reads

read.depth.tar=function(PlasmidDesigns){        # new function that requires 1 dataset (PlasmidDesigns)
  
  targetcount=c(PlasmidDesigns$match)           # vector with all targetcounts
  meantargets=mean(targetcount)                # compute mean for all targetcounts
  
  barplot(targetcount,xlab="targets",ylab="number of reads",main="read depth (targets)")     # plot number of reads vs. targets
  abline(h=meantargets,col="red",lty=2,lwd=2)                                                # add mean line
  legend("topright",c("mean"),cex=0.8,col="red",lty=2,lwd=2)                                 # and a legend
}

#quality assurance: plot number of genes vs. number of reads

read.depth.gene=function(PlasmidGenes){             # new function that requires 1 dataset (PlasmidDesigns)
  
    genecount=c(PlasmidGenes$match)                 # vector with all genecounts
  meangenes=mean(genecount)                         # compute mean for all genecounts
  
  barplot(genecount,xlab="genes",ylab="number of reads",main="read depth (genes)")       # plot number of reads vs. genes
  abline(h=meangenes,col="red",lty=2,lwd=2)                                              # add mean line
  legend("topleft",c("mean"),cex=0.8,col="red",lty=2,lwd=2)                              # and a legend
  
  gap.barplot(genecount,gap=c(21000,160000),xlab="genes",ylab="number of reads",main="read depth (genes)",ytics=c(0,5000,10000,15000,20000,165000,170000))    # plot number of reads vs. genes with gaped y-Axis
  abline(h=meangenes,col="red",lty=2,lwd=2)                                                                                                                       # add mean line
  legend("topleft",c("mean"),cex=0.8,col="red",lty=2,lwd=2)                                                                                                       # and a legend
}

#mean (per gene) vs. variance (per gene) 

mean.var.pergene=function(PlasmidDesigns,PlasmidGenes){    # new function that requires 2 datasets (PlasmidDesigns,PlasmidGenes)
  
  genecount=c(PlasmidGenes$match)                          # vector with all genecounts
  designspergene=c(PlasmidGenes$X..designs.total)          # vector with number of designs per gene
  meanpergene=c(NA,NA)
  
  for(i in 1:length(genecount)){                           # compute mean 
    meanpergene[i]=genecount[i]/designspergene[i]          # per design per gene
  }
                                                      
  designcount=c(PlasmidDesigns$match)                      # vector with all targetcounts
  varpergene=c(NA,NA)                                      #
  
  j=1                                                             #
  for(i in 1:length(genecount)){                                  # compute the variance
    varpergene[i]=var(designcount[j:(j+(designspergene[i]-1))])   # between all targets of one gene
    if(is.na(varpergene[i])){                                     #
      varpergene[i]=0                                             # set variance=0 if there is only one target
    }                                                             #
    j=(j+designspergene[i])                                       #
  }
  
  cormeanvar=cor(meanpergene,varpergene,method="pearson")         # compute corralation coefficient between mean and variance
  cormeanvar=round(cormeanvar,digits=4)                           # round
  
  plot(meanpergene,varpergene,xlab="mean (per gene)",ylab="variance (per gene)")    # plot variance per gene vs. mean per gene
  legend("topright",c("the correlation coefficient is:",cormeanvar),cex=0.6)        # add correlation coefficient to the plot
  #cat ("The pearson's correlation coefficient (mean per genes vs. variance per gene) is:",cormeanvar,"\n")
}

read.distr.targets(PlasmidDesigns)                #
read.distr.genes(PlasmidGenes)                    #
read.depth.tar(PlasmidDesigns)                    # execute all the functions
read.depth.gene(PlasmidGenes)                     #
mean.var.pergene(PlasmidDesigns,PlasmidGenes)     #
