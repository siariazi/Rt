# Rt
Files for generating data with OOPidemic package, estimate Rt from phylogeny with BEAST using BD skylihne plot or case count data using EpiEstim. \
EpiEstim_intro.R: An introductory script to run EpiEstim package (from the website). \
R0_v1.R: This is Chris’s OOpidemic code to run a simulation 1000 times and estimate R0. \
R0_v1_2.R: This is like R0_v1.R but no loop to see what’s going on. \
seir_v1.R: Chris’s code: SEIR model with OOpidemic generating fasta file and meta data. \
group_vingnett.rmd: Chris’s vingnett to show the package OOpidemic. \
seq_over_time.R: plotting the number of sequences over time from a fasta file. \
add_date_before.R: using metadata from OOPidemic to add date to fasta file. In this file I flip the dates from simulation, means that if time is 0 I convert it to max_time-0, when using fasta file generated from this script, in BEAUTi when setting tip dates, “before the present” (2nd) option should be used. \
add_date_since.R: like add_date_before.R, when using fasta file generated from this script, in BEAUTi when setting tip dates, “since some time in the past” (1st) option should be used. I prefer this version. \
bdsky_example_cond.rmd: estimating Re from beast results. The script is written by Louis du Plessis. In this version we are not estimating the time of origin and tree height should be given from opening the log file in tracer. \
bdsky_example_regular.rmd: estimating Re from beast results. The script is written by Louis du Plessis. In this version we are estimating the time of origin (I prefer this version). \
Workflow: metadata and fasta file from OOPidemic package should be combined together (add_data_since.R), then run the BESAT using the config explanined below and finally post-analysis using bdsky_example_regular.rmd. \ \
Configuring the BEAUti file: \
In setting Tip Dates be careful to choose since some time in the past \ 
Site model: subst model: HKY \
Go to priors: \
Birth Death Skyline Serial \
Initialization: \
Uncheck the estimate box for becomeUnifectousRate (I skip this unless I'm sure about the recovery rate) \
dimension of Reproductive number: 10 \
dimension of samplingProportion parameters to 10 \
Prior: \
Expand the options for Tree.t:alpha: Enter 36.5 (or whatever is recover_rate) for Become Uninfectious Rate (again I skip this) \
Set a Log Normal prior for both origin parameters with M = -0.5 and S = 0.2 \
reproductiveNumber_BDSKY: lognormal with M: 0.8 and S: 0.5 \
samplingProportion_BDSKY: beta distn. Alpha: 2 beta: 1000 












