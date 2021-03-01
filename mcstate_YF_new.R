library(mcstate)
sero_data=read.csv("sero_data.csv",header=TRUE) #Seroprevalence data for comparison
source("compare_function_YF.R")
sero_data2=data.frame(year=c(1:length(sero_data$year)),positives=sero_data$positives)
sero_data3 <- particle_filter_data(data = sero_data2,time = "year",rate = 1)
#-----------------------------------------------------------------------------------------------------
model_file = "YFOutbreak_dust.cpp"
YFmodel <- dust::dust(model_file) #Compile model
#-----------------------------------------------------------------------------------------------------
n_particles <- 4 #If this is set to more than 1, pmcmc() returns error
n_threads <- 4 #This can be set to more than 1 without causing an error
filter <- particle_filter$new(data = sero_data3,model = YFmodel,
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
n_steps=1000
control <- pmcmc_control(n_steps=n_steps,save_state = TRUE,save_trajectories = TRUE,progress = TRUE)
#-----------------------------------------------------------------------------------------------------
pmcmc_run <- pmcmc(mcmc_pars, filter, control = control)

probs=pmcmc_run$probabilities[,3]
plot(probs,xlab="Step",ylab="log posterior likelihood")

steps_select=c(0,n_steps/4,n_steps/2,n_steps)+1
age_min=0
age_max=18
xvalues=pmcmc_run$trajectories$step
obs_values=sero_data$positives/sero_data$samples
sero_out1=sero_out2=sero_out3=sero_out4=rep(0,length(xvalues))
for(i in 1:length(xvalues)){
  sero_out1[i]=mean(pmcmc_run$trajectories$state[c((age_min+1):age_max)+1,steps_select[1],i])
  sero_out2[i]=mean(pmcmc_run$trajectories$state[c((age_min+1):age_max)+1,steps_select[2],i])
  sero_out3[i]=mean(pmcmc_run$trajectories$state[c((age_min+1):age_max)+1,steps_select[3],i])
  sero_out4[i]=mean(pmcmc_run$trajectories$state[c((age_min+1):age_max)+1,steps_select[4],i])
}
ylabel=paste("Seroprevalence ages ",age_min,"-",age_max,sep="")
ymax=max(sero_out1,sero_out2,sero_out3,sero_out4,obs_values)

matplot(x=xvalues,y=sero_out1,xlab="Year",ylab=ylabel,type="l",lty=1,col=2,ylim=c(0,ymax))
matplot(x=xvalues,y=sero_out2,type="l",col=3,lty=2,add=TRUE)
matplot(x=xvalues,y=sero_out3,type="l",col=4,lty=3,add=TRUE)
matplot(x=xvalues,y=sero_out4,type="l",col=5,lty=4,add=TRUE)
matplot(x=sero_data$year-2000,obs_values,type="b",pch=1,col=1,add=TRUE)
legend("topright",legend=c(paste("Step ",steps_select[1],sep=""),paste("Step ",steps_select[2],sep=""),
                           paste("Step ",steps_select[3],sep=""),paste("Step ",steps_select[4],sep=""),
                           "Observed"),col=c(2,3,4,5,1),lty=c(1,2,3,4,1))

FOI_final=exp(pmcmc_run$pars[n_steps+1,1])[[1]]
R0_final=exp(pmcmc_run$pars[n_steps+1,2])[[1]]
FOI_final
R0_final