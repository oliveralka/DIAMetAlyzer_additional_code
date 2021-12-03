concentrationseries.Rmd:
Rmd script for the processing of quantitative data based the analyis of the pesticide mix using the 
DIAMetAlyzer pipeline. The script used to generate the figures representing the quantification results. 

Confusion_metrics.Rmd:
Rmd script to plot the graphics based on the multi confusion matrix (see comparison).
Here, FDR (ground truth) vs estimated FDR, ROC, as well es precision, recall plots have been
created. 

lod.Rmd:
Rmd script calculates the mean intensity over the dilution series and 
the signal-to-noise based on the last detected signal (regarded as noise).
With this, we can assess the limit of detection (LOD). 