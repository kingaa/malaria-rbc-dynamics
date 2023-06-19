# Part I: Final files for fitting a hierarchical model to NW11 data using penalised regression 

## Model types (x)
- Full model: fits all five boxes (including control mice) and does not attempt to fit missing data
- 3 boxes: fits the control box, pABA = 0% and pABA = 0.5% boxes and does not attempt to fit missing data
- 4 boxes: fits all boxes except pABA = 5% and DOES attempt to fit missing data

## File types
- FDA_x.R: R script to run fitting for model specified by x
- FDA_x_results.Rmd: R Markdown file to import and visualise model fitting results; also reads in FDA_x.rda to give information on convergence and final likelihood
- FDA_x_results.html: HTML output of FDA_x_results.Rmd
- parameters_results_x.csv: parameter estimates for model specified by x
- trajectory_results_x.csv: trajectory estimates for model specified by x
- trajectory_guesses.csv: provides trajectory guesses for three boxes (pABA = 0%, 0.05% and 0.5%); superseded by using data and model to obtain initial trajectory guesses

# Part II: Files for obtaining group level trajectories from Wale et al. 2017 results
## File types
- POMP_GroupLevel.R: initial R script for obtaining group-level trajectories (superseded by POMP_GroupLevel.Rmd)
- POMP_GroupLevel.Rmd: R Markdown for obtaining and visualising group-level trajectories
- POMP_GroupLevel.html: HTML output of POMP_GroupLevel.Rmd
- m5pf4.rds: used to obtain MLEs for sigmas, initial values, betas and doses
- m5sm1.rds: individual replicates/trajectories read into R and processed to obtain group-level trajectories
