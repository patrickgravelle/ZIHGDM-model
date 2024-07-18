# ZIHGDM-model

The following files allow for the construction of new ZIHGDM models. There are three sets of files, those completed for TDP-43, G85R, and 5XFAD data. None of the data is included here.

Constructing/defining the model specifications are done jointly through the ZIHGDM_model_Stan_file_XX.stan and ZIHGDM_model_Runfile_XX.R files. The output of the ZIHGDM_model_Runfile_XX.R file can be used to generate posterior predictive plots with the posterior_predictive_plots_XX.R files to determine model fit to the data. The posterior means stratified by hour and behavior can then be generated from the posterior_mean_analysis_TDP43.R file. Each file should not be run without accounting for the corresponding changes in the model specifications.
