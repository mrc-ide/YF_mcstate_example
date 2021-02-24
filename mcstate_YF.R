library(mcstate)
setwd("C:/Users/Work_KJF82/Documents/0 - R files/Tests/MCState testing/mcstate YF example") #Set to folder containing files
model_file = "YFOutbreak_dust.cpp"
sero_data=read.csv("sero_data.csv",header=TRUE) #Seroprevalence data for comparison
source("compare_function_YF.R")

#Parameters currently set manually in CPP file
{
# pop_data=read.csv("population_data.csv",header=TRUE)
# F_immune_initial=0.8
# vaxrate_infant=F_immune_initial/365
# n_years=length(pop_data$year)
# n_age_brackets=ncol(pop_data)-1
# pop_data_vector=c()
# for(i in 1:n_age_brackets){
#   pop_data_vector=append(pop_data_vector,pop_data[[i+1]])
# }
}
#-----------------------------------------------------------------------------------------------------
YFmodel <- dust::dust(model_file) #Compile model
sero_data2 <- particle_filter_data(data = sero_data,time = "year",rate = 1) #Format seroprevalence data for comparison
n_particles <- 1 #If this is set to more than 1, pmcmc() returns error
n_threads <- 1 #This can be set to more than 1 without causing an error
filter <- particle_filter$new(data = sero_data2,model = YFmodel,
                              n_particles = n_particles, n_threads=n_threads, compare = sero_compare,seed = 1L)
#-----------------------------------------------------------------------------------------------------
# Parameters for fitting - log reproduction number and log external FOI
logR0 <- pmcmc_parameter("logR0", log(5.0), min = log(0.5),max=log(10.0),
                         prior = function(p)
                           dnorm(p, mean = log(5.0), sd=log(2.0), log = FALSE))
logfoi <- pmcmc_parameter("logfoi", log(1.0e-4), min = log(5.0e-6),max=log(1.0e-3),
                          prior = function(p)
                            dnorm(p, mean = log(1.0e-6), sd=log(5.0), log = FALSE))
proposal_matrix <- matrix(c(0.001, 0, 0, 0.001), nrow = 2, ncol = 2, byrow = TRUE)
pars=list(logfoi=logfoi,logR0=logR0)
mcmc_pars <- pmcmc_parameters$new(pars, proposal_matrix)
#-----------------------------------------------------------------------------------------------------
control <- pmcmc_control(n_steps=100,save_state = TRUE,save_trajectories = TRUE,progress = TRUE)
pmcmc_run <- pmcmc(mcmc_pars, filter, control = control)

probs=pmcmc_run$probabilities[,3]
plot(probs,xlab="Step",ylab="log posterior likelihood")
