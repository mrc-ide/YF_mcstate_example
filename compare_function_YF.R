sero_compare <- function(state, observed, pars = NULL) {
  
  n_lines=length(observed$positives)
  sero_modelled <- state[, , drop = TRUE]
  sero_modelled2=rep(0.5,n_lines)

  for(i in 1:n_lines){
    sero_modelled2[i]=mean(sero_modelled[c(1:18)])
  }

  LogLikelihood = sum(lgamma(observed$samples+1)-lgamma(observed$positives+1)
                      -lgamma(observed$samples-observed$positives+1) +
                        observed$positives*log(sero_modelled2) +
                        (observed$samples-observed$positives)*log(1.0-sero_modelled2))
  
  #LogLikelihood=runif(1,-10,-1)
  return(LogLikelihood)
  
}