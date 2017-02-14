n=commandArgs(trailing=T)

doSingleMarkerTests(data)
{
}

doPerms(data,nperms)
{
}

simSites=function(nsites,nsam,a,b)
    #' Site frequencies are modeled as beta-distributed
    #' with coefficients a and b
{
    return(floor(rbeta(nsites,a,b)*nsam)+1)
}


nsam=n[1]
nsites=n[2]
nperms=n[3]
nreps=n[4]
beta_a=n[5]
beta_b=n[6]

