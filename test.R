source("sim.R")

g=simSites(1000,500,1,5)

for (noise in c(0.1,0.25,0.5,1))
{
    p=simTraits(g,noise)
    pv=doSingleMarkerTests(g,p)
    obs=getMedianDiffs(g,p,1:ncol(p))
    bounds=doPerms(g,p,100)
    print(noise)
    #print(bounds)
    #print(paste(min(obs),max(obs)))
    #print(is.numeric(bounds$lower))
    #print(paste(length(which(obs <= bounds$lower)),length(which(obs>=bounds$upper))))
    sig_perm_low = which(obs <= bounds$lower)
    sig_perm_hi = which(obs >= bounds$upper)
    pv_low = pv[sig_perm_low]
    pv_hi = pv[sig_perm_hi]
    #Print true positive rate for permutation method
    print(paste(length(which(pv_low<=0.05)),length(pv_low)))
    print(paste(length(which(pv_hi<=0.05)),length(pv_hi)))
    #Print what comes out of single-marker testing
    print(length(which(pv<=0.05)))
    #Overlap b/w single-marker p-values and what perms find
    print(length(intersect(which(pv<=0.05),c(sig_perm_hi,sig_perm_low))))
}
