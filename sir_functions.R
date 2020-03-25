## As above, but each age group has its own transmission parameter
general_sir <- function(t,y, pars, C, Nage,Nimmunity){
  beta <- pars[1] ## Overall transmission rate
  Tg <- pars[length(pars)] ## Recovery time
  alphas <- pars[2:(length(pars)-1)] ## Vector of susceptibility parameters
  
  ## Generate matrix for compartment sizes
  sir <- matrix(y,ncol=Nage*Nimmunity,nrow=3)
  dS <- -(alphas*sir[1,]) * (beta*sir[2,] %*% t(C))
  dR <- sir[2,]/Tg
  dI <- -dS - dR
  
  tmp <- as.vector(rbind(dS,dI,dR))
  
  return(list(c(tmp)))
}

## As above, but each age group has its own transmission parameter
general_sir_explicit <- function(t,y, pars, C, Nage,Nimmunity){
  beta <- pars[1] ## Overall transmission rate
  Tg <- pars[length(pars)] ## Recovery time
  alphas <- pars[2:(length(pars)-1)] ## Vector of susceptibility parameters
  
  ## Generate matrix for compartment sizes
  sir <- matrix(y,ncol=Nage*Nimmunity,nrow=3)
  N <- colSums(sir)
  M <- t(apply(C, 1, function(x) x/N))
  M[!is.finite(M)] <- 0
  dS <- -(alphas*sir[1,]) * apply(M, 1, function(x) sum(beta*sir[2,]*x))
  dR <- sir[2,]/Tg
  dI <- -dS - dR
  
  tmp <- as.vector(rbind(dS,dI,dR))
  
  return(list(c(tmp)))
}

#' Epidemic Final Size Calculation ODE
#' 
#' Calculates the final size of an epidemic given 2-dimensional population categorisation eg. age and immunity class using an SIR model
#' @param C1 the normalised contact matrix of contact frequencies between each age and immunity class
#' @param beta the disease specific beta (transmission rate). Note that another parameter will mediate the contact rate
#' @param Tg number of days spent infectious (ie. 1/gamma)
#' @param Ns the matrix of population sizes for each age/immunity combination (non-normalised) (ie. rows = ages, cols = immunity classes)
#' @param alphas a vector of values between 0 and 1 matching the number of immunity classes
#' @param ts vector of times to solve over
#' @param beta_scales vector giving relative transmissibility (0-1) for each age class
#' @param age_seed which age group to seed in
#' @param immunity_seed which immunity class to seed in
#' @return an NxM matrix of attack rates (ie. proportion of susceptibles becoming infected)
#' @seealso \code{\link{epi_final_size}}
#' @export
epi_ode_size <- function(C1, beta, Tg, Ns, alphas, 
                                ts=seq(1,365,by=1),
                                age_seeds=c(1),immunity_seed=1,
                                ver="fast",return_peak=FALSE){
  C <- C1
  Nage <- nrow(Ns)
  Nimmunity <- ncol(Ns)
  
  ## Generate starting populations.
  ## This generates a vector with S, I and R for each age class and immunity class
  ## ie. S_11, I_11, R_11, S_12, I_12, R_12, S_22, ... etc.
  ## where S_ak gives the number susceptible in age class 1, immunity class k
  long_Ns <- as.numeric(t(Ns))
  start <- NULL
  start[1] <- long_Ns[1]
  start[2] <- 0
  start[3] <- 0
  index <- 4
  for(i in 2:length(long_Ns)){
    start[index] <- long_Ns[i]
    index <- index + 1
    start[index] <- 0
    index <- index + 1
    start[index] <- 0
    index <- index + 1
  }
  start[start < 0] <- 0
  
  ## Move one susceptible to the infected classes for the seed population
  for(age_seed in age_seeds){
    start[(age_seed-1)*3*Nimmunity + (immunity_seed-1)*3+1] <- start[(age_seed-1)*3*Nimmunity + (immunity_seed-1)*3+1] - 1
    start[(age_seed-1)*3*Nimmunity + (immunity_seed-1)*3+2] <- start[(age_seed-1)*3*Nimmunity + (immunity_seed-1)*3+2] + 1
  }
  #start[(age_seed + (immunity_seed-1) - 1)*3 + 1] <- start[(age_seed + (immunity_seed-1) - 1)*3 + 1] - 1
  #  start[(age_seed + (immunity_seed-1) - 1)*3 + 2] <- start[(age_seed + (immunity_seed-1) - 1)*3 + 2] + 1
  
  #beta_scales <- rep(beta_scales, each=Nimmunity)
  if(ver == "fast"){
    y <- ode(y=start,t=ts,func=general_sir, parms=c(beta,alphas,Tg),C=C,Nage=Nage,Nimmunity=Nimmunity)
  } else {
    y <- ode(y=start,t=ts,func=general_sir_explicit, parms=c(beta,alphas,Tg),C=C,Nage=Nage,Nimmunity=Nimmunity)
  }
  ## Pull out the solved model.
  y <- as.data.frame(y)
  y <- y[,2:ncol(y)]
  
  ## Label which disease state each column is
  colnames(y) <- rep(c("S","I","R"), Nage*Nimmunity)
  
  ## Pull out the recovered population to get final size
  recovered <- y[,which(colnames(y) == "R")]
  prevalence <- y[,which(colnames(y) == "I")]
  total_prev <- rowSums(prevalence)
  peak_time <- ts[which.max(total_prev)]
  ## Final size is number in recovered at end minus number recovered at start
  final_incidence <-  recovered[nrow(recovered),] - recovered[1,]
  final_incidence <- as.numeric(final_incidence/long_Ns)
  final_incidence[is.nan(final_incidence)] <- 0
  A <- matrix(final_incidence,nrow=Nage,ncol=Nimmunity,byrow=T)
  if(!return_peak){
    return(A)
  } else {
    return(list(A, peak_time))
  }
}
