# plot 2 datasets against each other

generead1.generead2=function(SelectionGenes1,SelectionGenes2){       # new function that requires 2 datasets (SelectionGenes1,SelectionGenes2)
  
  genecount1=c(SelectionGenes1$match)                                # vectors with genecounts
  genecount2=c(SelectionGenes2$match)                                # for both datasets
  
  plot(genecount1,genecount2,pch=20,xlab="reads per gene (dataset 1)",ylab="reads per gene (dataset 2)")    # plot genecounts of both datasets against each other
  
  genecount1[length(genecount1)]=0                      # "delete" the last datapoint which
  genecount2[length(genecount2)]=0                      # lies far away from all others
  
  corgeneread=cor(genecount1,genecount2)                # compute correlation coefficient            
  corgeneread=round(corgeneread,digits=4)               # round
  
  plot(genecount1,genecount2,pch=20,xlab="reads per gene (dataset 1)",ylab="reads per gene (dataset 2)")    # plot genecounts of both datasets against each other
  legend("topleft",c("the correlation coefficient is:",corgeneread),cex=0.75)                               # add correlation coefficient to plot
}

meandens1.meandens2=function(SelectionGenes1,SelectionGenes2){       # new function that requires 2 datasets (SelectionGenes1,SelectionGenes2)
                                  
  genecount1=c(SelectionGenes1$match)                                # vector with all genecounts
  designspergene1=c(SelectionGenes1$X..designs.total)                # vector with number of designs per gene
  meanpergene1=c(NA,NA)                                              #
  
  for(i in 1:length(genecount1)){                                    # compute mean per design per gene
    meanpergene1[i]=genecount1[i]/designspergene1[i]                 # for dataset 1
  }
  
  meandensity1=density(meanpergene1)                                 # compute density function of mean per gene
  
  plot(meandensity1,ylim=c(0,0.0030),col="red",xlab="mean per gene",ylab="density",main="density distribution")   # plot density function of mean per gene for dataset 1
                                                              
  genecount2=c(SelectionGenes2$match)                                # vector with all genecounts
  designspergene2=c(SelectionGenes2$X..designs.total)                # vector with number of designs per gene
  meanpergene2=c(NA,NA)                                              # 
 
  for(i in 1:length(genecount2)){                                    # compute mean per design per gene
    meanpergene2[i]=genecount2[i]/designspergene2[i]                 # for dataset 2
  }
  
  meandensity2=density(meanpergene2)                                 # compute density function of mean per gene
  
  lines(meandensity2,col="blue")                                                  # add density function of mean per gene for dataset 2 to plot of data 1
  legend("topright",c("run 1","run 2"),col=c("red","blue"),lty=1,lwd=1,cex=0.75)  # add a legend
}

generead1.generead2(SelectionGenes1,SelectionGenes2)     # execute all the functions
meandens1.meandens2(SelectionGenes1,SelectionGenes2)     #
