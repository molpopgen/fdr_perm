library(qvalue) #from bioconductor

doSingleMarkerTests=function(genotypes,phenotypes)
    #' This is the correct p-value, which is the regression
    #' of standardized phenotype onto genotype, on a per-variant basis
{
    rv=array()
    for(i in 1:ncol(phenotypes))
    {
        #Turn phenotype into z-score
        np=(phenotypes[,i]-mean(phenotypes[,i]))/sd(phenotypes[,i])
        rv[i]=summary(glm(np~genotypes[,i]))$coefficients[,4][2]
    }
    return(rv)
}

getMedianDiffs=function(genotypes,phenotypes,labels)
{
    rv=array()
    for(i in 1:ncol(phenotypes))
    {
        without=median(phenotypes[which(genotypes[,i][labels]==0),i])
        with=median(phenotypes[which(genotypes[,i][labels]==1),i])
        rv[i]=with-without
    }
    return(rv)
}

doPerms=function(genotypes,phenotypes,nperms)
{
    labels=1:nrow(genotypes)
    rv=NA
    for(i in 1:nperms)
    {
        labels=sample(labels)
        md=getMedianDiffs(genotypes,phenotypes,labels)
        if(i==1)
        {
            rv=md
        }
        else
        {
            rv=c(rv,md)
        }
    }
    return(list(lower=as.numeric(quantile(rv,0.025)),upper=as.numeric(quantile(rv,0.975))))
}

simTraits=function(data,noise,tpr,upreg)
    #' Mean expression level is 1 + noise for absence 
    #' genotype and upreg + noise for presence genotype
    #' Noise is a Gaussian
{
    phenotypes=data
    truePositives=array()
    for(i in 1:ncol(phenotypes))
    {
        phenotypes[,i][which(data[,i]==0)] = 1. 
        if(runif(1)<=tpr)
        {
        phenotypes[,i][which(data[,i]==1)] = upreg 
        truePositives[i]=1
        }
        else {
        phenotypes[,i][which(data[,i]==1)] = 1. 
        truePositives[i]=0
        }
        #Add the noise term for each individual.
        #This reflects cis/trans/E stuff
        for(j in 1:nrow(phenotypes))
        {
            phenotypes[j,i]=phenotypes[j,i]+rnorm(1,mean=0,sd=noise)
        }
    }
    return(list(phenotypes=phenotypes,truePositives=truePositives))
}

simSites=function(nsites,nsam,a,b)
    #' Site frequencies are modeled as i.i.d. & beta-distributed
    #' with coefficients a and b
{
    sampleCounts=floor(rbeta(nsites,a,b)*nsam)+1
    data=matrix(nrow=nsam,ncol=nsites,data=0.0)
    for(i in 1:length(sampleCounts))
    {
        data[,i][sample(nsam,sampleCounts[i],FALSE)]=1
    }
    return(data)
}

doStudy=function(nsites,nsam,noise,tpr,upreg,a,b)
{
    g=simSites(nsites,nsam,a,b)
    ptp = simTraits(g,noise,tpr,upreg)
    p=ptp$phenotypes
    truePositives=which(ptp$truePositives==1)
    pv=doSingleMarkerTests(g,p)
    qv=qvalue(pv,fdr=0.05)  #Storey's FDR method at 0.05
    qvsig=which(qv$sigificant == TRUE)
    obs=getMedianDiffs(g,p,1:ncol(p))
    bounds=doPerms(g,p,100)
    sig_perm_low = which(obs <= bounds$lower)
    sig_perm_hi = which(obs >= bounds$upper)
   
    fdr_qv = length(setdiff(qvsig,truePositives))/length(which(pv<=0.05))
    fdr_sm = length(setdiff(which(pv<=0.05),truePositives))/length(which(pv<=0.05))
    tpr_sm = length(intersect(which(pv<=0.05),truePositives))/length(which(pv<=0.05))
    fdr_perm = length(setdiff(c(sig_perm_hi,sig_perm_low),truePositives))/(length(sig_perm_hi)+length(sig_perm_low))
    tpr_perm = length(intersect(c(sig_perm_hi,sig_perm_low),truePositives))/(length(sig_perm_hi)+length(sig_perm_low))
    return(list(fdr_sm=fdr_sm,fdr_sm_qv=fdr_qv,fdr_perm=fdr_perm,tpr_sm=tpr_sm,tpr_perm=tpr_perm))
}
#n=commandArgs(trailing=T)
#nsam=n[1]
#nsites=n[2]
#nperms=n[3]
#nreps=n[4]
#beta_a=n[5]
#beta_b=n[6]

