# smc.covid
smc.covid is an R package related to the paper A sequential Monte Carlo approach to estimate a time varying reproduction number in infectious disease models: the COVID-19 case by Storvik et al (2021).

The package can be installed in R by the commands

*library(devtools)*

*install.packages("spread", repos =  "https://folkehelseinstituttet.github.io/drat/")*

*install_git("https://github.com/geirstorvik/smc.covid")*


The main algorithm is a version of sequential Monte Carlo (the simplest bootstrap filter) where there are two underlying processes, one defining the daily reproduction numbers and one defining a SEIR model.

The main routine is **smc_loglik_mem** which runs the sequential Monte Carlo algorithm. Results are partly given in a returning object and partly written to files in a specified result directory.

There are also two routines for performing particle PMCMC. **smc_PMCMC_par** updates the a and sigma parameters in the AR_mu_0 model while **smc_PMCMC** is somewhat more general (but currently not working!)

An example dataset based on real data from Norway, **Norway_adj**, is provided together with a script, **run_Norway.R** (available in the test directory, you probably need to download this file separatly). Note that in order to obtain reasonable results at least 5000 particles should be used, but a run is time-consuming. The implementation is parallelized, so using 10-50 cores is preferable. For testing, smaller number of particles and/or earlier end-dates can be used.

For testing of installation and that it works, an alternative script, **run_Norway_test.R** is also given which uses only 4 cores, have much less particles and only run for a short time period. Results should not be trusted here though.





