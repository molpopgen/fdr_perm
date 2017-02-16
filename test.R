library(tibble)
library(dplyr)
source("sim.R")
NSIMS_PER_PARAM=100

sim_results=tibble(upreg=numeric(),tpr=numeric(),noise=numeric(),fdr_sm=numeric(),fdr_sm_qv=numeric(),fdr_perm=numeric(),tpr_sm=numeric(),tpr_perm=numeric())
for(upreg in c(2,5,10))
    {
    for (tpr in seq(0.05,0.2,0.05))
        {
        for (noise in c(0.1,0.25,0.5,1))
        {
            for(replicate in 1:NSIMS_PER_PARAM)
            {
                results = doStudy(1000,20,noise,tpr,upreg,1,5)
                sim_results=add_row(sim_results,upreg,tpr,noise,fdr_sm=results$fdr_sm,fdr_sm_qv=results$fdr_sm_qv,fdr_perm=results$fdr_perm,tpr_sm=results$tpr_sm,tpr_perm=results$tpr_perm)
            }
        }
    }
}

means = sim_results %>% 
    group_by(upreg,tpr,noise) %>%
    summarise(mfdr_sm=mean(fdr_sm),mfdr_sm_qv=mean(fdr_sm_qv),mfdr_perm=mean(fdr_perm))#,mtpr_sm=mean(tpr_sm),mtpr_perm=mean(tpr_perm))

write.table(data.frame(means),file="mean_fdr.txt",quote=F,row.names=F)
