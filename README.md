# Github repo for paper "Improving Power  through Adjustment for Prognostic Variables in Group Sequential Trial Designs: Impact of Baseline Variables, Short-Term Outcomes, and Treatment Effect Heterogeneity"

(a.k.a. Group Sequential Design with TMLE)

Authors: Tianchen Qian, Michael Rosenblum, Huitong Qiu

Link to working paper: https://biostats.bepress.com/jhubiostat/paper285/



## To reproduce simulation results in the paper:

### 1. "Run main_simtrial.R"

This file is currently set up so that one needs to submit a batch array job (with job id 1-40) to a cluster. The main_simtrial.slurm is what I used to submit batch job. One can easily modify the code in main_simtrial.R to run the entire simulation sequentially, and it would take about 160 hours to complete on a single core 3.0 GHz CPU.

The simulation depends on a data set that consists of 100 patients data from MISTIE Phase II trial. Unfortunately, the data cannot be made public.

### 2. Run "main_make_table.R"

This file collects output from the above parallelized main_simtrial.R, and computes quantities for the following tables in the paper: Table 3, Table 4, Table D.2, Table G.2, Table G.3.

### 3. Run "make plot - GST_TMLE paper.R"

This file makes the figures in the paper: Figure 1, Figure 2, Figure F.1, Figure F .2.

## To reproduce "ARE from theory" in Table 3, Table G.1, Table G.2:

Run main_calculate_RE_from_dgm.R.


## Structure of the code

* auxiliary_code_construct_generative_model/: auxiliary code that was used to find tuning parameters for the generative model (a.k.a. DGM, data generating mechanism)

* auxiliary_code_determine_nmax/: auxiliary code that was used to find the maximum sample size (nmax) to achieve 80% power for the simulated trials under various prognostic settings.

* parallel_seed/: random seeds for parallel jobs.

* source/: function definition used in the simulation

  * fcn_DGM.R: generative model
  * fcn_compcov.R: calculate (empirical) covariance matrix of the test statistics from all simulated trials
  * fcn_comppower.R: calculate type I error, power, etc. for a simulated trial
  * fcn_err_sp_bdry.R: calculate the testing boundaries using error spending approach
  * fcn_evaltrial.R: calculate the percentage of early stopping at each stage
  * fcn_reverse_prob.R: calculate the probability of "reverse" (favoring H1 at interim analysis but fail to reject at decision analysis, or vice versa)
  * fcn_simtrials.R: simulate trials to get test statistics
  * fcn_trialresults.R: some utility functions used in other functions
  * v.R: configuration of a trial, including error spending functions, etc.
