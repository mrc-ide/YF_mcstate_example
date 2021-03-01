sero_compare <- function(state, observed, pars = NULL) {
  exp_noise <- 1e6
  prevalence_modelled <- state[6, , drop = TRUE]
  for(i in 7:23){
    prevalence_modelled=prevalence_modelled+state[7,,drop=TRUE]
  }
  prevalence_modelled=prevalence_modelled/18
  prevalence_modelled=as.integer(prevalence_modelled*100000)
  
  prevalence_observed <- observed$positives*100
  #prevalence_modelled <- prevalence_observed+sample(c(1:5))[1]
  lambda <- prevalence_modelled +
    rexp(n = length(prevalence_modelled), rate = exp_noise)
  dpois(x = prevalence_observed, lambda = lambda, log = TRUE)
}


# sero_compare <- function(state, observed, pars = NULL) {
#   
#   n_lines=length(observed$positives)
#   sero_modelled <- state[, , drop = TRUE]
#   sero_modelled2=rep(0.5,n_lines)
# 
#   for(i in 1:n_lines){
#     sero_modelled2[i]=mean(sero_modelled[c(2:19)])
#   }
# 
#   LogLikelihood = sum(lgamma(observed$samples+1)-lgamma(observed$positives+1)
#                       -lgamma(observed$samples-observed$positives+1) +
#                         observed$positives*log(sero_modelled2) +
#                         (observed$samples-observed$positives)*log(1.0-sero_modelled2))
#   
#   #LogLikelihood=runif(1,-10,-1)
#   return(LogLikelihood)
#   
# }