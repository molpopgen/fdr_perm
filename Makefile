all: mean_fdr.txt fdr.pdf

mean_fdr.txt: sim.R test.R
	R CMD BATCH test.R

fdr.pdf: mean_fdr.txt plot.R
	R CMD BATCH plot.R
