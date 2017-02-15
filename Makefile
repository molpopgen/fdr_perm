all: mean_fdr.txt

mean_fdr.txt: sim.R test.R
	R CMD BATCH test.R
