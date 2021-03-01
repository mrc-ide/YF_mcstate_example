source("compare_function_YF.R")
sero_data2=data.frame(year=c(1:length(sero_data$year)),positives=sero_data$positives)
sero_data3 <- particle_filter_data(data = sero_data2,time = "year",rate = 1)

sero_compare(state=pmcmc_run$state, observed=sero_data3, pars = NULL)


exp_noise <- 1e6
prevalence_modelled <- pmcmc_run$state[2, , drop = TRUE]
for(i in 3:19){
  prevalence_modelled=prevalence_modelled+pmcmc_run$state[i,,drop=TRUE]
}
prevalence_modelled=prevalence_modelled/18
prevalence_modelled=as.integer(prevalence_modelled*1000)

prevalence_observed <- observed$positives
#prevalence_modelled <- prevalence_observed+sample(c(1:5))[1]
lambda <- prevalence_modelled + rexp(n = length(prevalence_modelled), rate = exp_noise)
dpois(x = prevalence_observed, lambda = lambda, log = TRUE)