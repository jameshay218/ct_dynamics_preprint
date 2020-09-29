# ct_inference_preprint
 Code required to regenerate the figures and analyses in Hay & Kennedy-Shaffer et al.

## General framework
- Data are stored in the data folder.
- R functions are stored in the code folder.
- All code to run analyses and generate figures is in the scripts folder.
- Model parameters are typically controlled by the corresponding table in the pars folder. Typically, the file with "fitted" in the name is after the optim fitting procedure. Change values, MCMC control and bounds in these tables.
- Most normal figures get saved to the figs folder.
- Figures generated as part of the fitting procedure tend to get saved to the results folder, in a sub-folder with name matching the analysis.

## lazymcmc
The lazymcmc package is needed throughout this project. It can be installed with `devtools::install_github("jameshay218/lazymcmc")`

## Working directory
The first thing that needs changing in each script is the starting working directory at the top of each script.

## Priors for fitting
Check `code/priors.R`. The priors for model fitting are currently hard-coded.

## fit_viral_load_model.R
- Uses optim to find the viral kinetics parameters that give a reasonable least squares. 
- It fits the old gamma model, the original single hinge function and the new double-hinge function. 
- Note that optim is quite sensitive to starting conditions, so find the lines that call optim and check that the starting parameters give fairly reasonable trajectories.