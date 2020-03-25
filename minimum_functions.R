



epi_ode_size_base <- function(C1, beta, Tg, Ns, ts=seq(1,365,by=1), prop_immune){
  C <- C1
  Nage <- nrow(Ns)
  Ntiter <- ncol(Ns)
  
  long_Ns <- as.numeric(t(Ns))
  start <- NULL
  start[1] <- (1-prop_immune[1])*long_Ns[1]-1
  start[2] <- 1
  start[3] <- prop_immune[1]*long_Ns[1]
  index <- 4
  for(i in 2:length(long_Ns)){
    start[index] <- (1-prop_immune[i])*long_Ns[i]-1
    index <- index + 1
    start[index] <- 1
    index <- index + 1
    start[index] <- long_Ns[i]*prop_immune[i]
    index <- index + 1
  }
  start[start < 0] <- 0
  
  y <- ode(y=start,t=ts,func=general_sir_base, parms=c(beta,alphas,Tg),C=C,Nage=Nage,Ntiter=Ntiter)
  y <- as.data.frame(y)
  y <- y[,2:ncol(y)]
  colnames(y) <- rep(c("S","I","R"), Nage*Ntiter)
  #return(y)
  recovered <- y[,which(colnames(y) == "R")]
  
  final_incidence <-  recovered[nrow(recovered),] - recovered[1,]
  final_incidence <- as.numeric(final_incidence/long_Ns)
  
  #A <- NULL
  #for(i in 1:((ncol(y)-1)/3)){
  #  A[i] <- y[nrow(y),(i-1)*3+4]/y[1,(i-1)*3+2]
  #}
  A <- matrix(final_incidence,nrow=Nage,ncol=Ntiter,byrow=T)
  
  return(A)
}

#' Epidemic Final Size Calculation ODE
#' 
#' Calculates the final size of an epidemic given 2-dimensional population categorisation eg. age and immunity class using an SIR model
#' @param C1 the non-normalised contact matrix of contact frequencies between each age class
#' @param beta the disease specific beta (transmission rate). Note that another parameter will mediate the contact rate
#' @param Tg number of days spent infectious (ie. 1/gamma)
#' @param Ns the matrix of population sizes for each age/titre combination (non-normalised) (ie. rows = ages, cols = immunity classes)
#' @param alphas a vector of values between 0 and 1 matching the number of immunity classes
#' @return an NxM matrix of attack rates (ie. proportion of susceptibles becoming infected)
#' @seealso \code{\link{epi_final_size}}
#' @export
epi_ode_size <- function(C1, beta, Tg, Ns, alphas, ts=seq(1,365,by=1)){
  #C <- setup_C(C1, Ns)
  C <- C1
  Nage <- nrow(Ns)
  Ntiter <- ncol(Ns)
  
  long_Ns <- as.numeric(t(Ns))
  start <- NULL
  start[1] <- long_Ns[1]-1
  start[2] <- 1
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
  
  y <- ode(y=start,t=ts,func=general_sir, parms=c(beta,alphas,Tg),C=C,Nage=Nage,Ntiter=Ntiter)
  y <- as.data.frame(y)
  y <- y[,2:ncol(y)]
  colnames(y) <- rep(c("S","I","R"), Nage*Ntiter)
  recovered <- y[,which(colnames(y) == "R")]
  
  final_incidence <-  as.numeric(recovered[nrow(recovered),] - recovered[1,])
  final_incidence <- as.numeric(final_incidence/long_Ns)
  
  #A <- NULL
  #for(i in 1:((ncol(y)-1)/3)){
  #  A[i] <- y[nrow(y),(i-1)*3+4]/y[1,(i-1)*3+2]
  #}
  A <- matrix(final_incidence,nrow=Nage,ncol=Ntiter,byrow=T)
  
  return(A)
}



# Make a parametric smooth WAIFW
# from Farrington, J. Amer. Stat. Assoc.; 2005, 100 p370;
#
#  parameters -
#     age class boundries - the upper age limit for each age class in years,
#     parameters defining the shape:
#        mu: age of highest contact increases with mu
#        gam: width around equal age diagonal increases with gam
#        sig: decreases strength in other diagonal (shrinks high trans int)
#        delta:  background homogeneous contact rate
#
#Returns -
#   a smooth WAIFW matrix

get.smooth.WAIFW<-function(age.class.boundries = (1:120/12),
                           mu=12.71,sig=0.69, gam=0.17, delta=0){
  
  n.age.cats <- length(age.class.boundries)
  ages.to.use <- (age.class.boundries +
                    c(0,age.class.boundries[2:n.age.cats-1]))/2
  
  #gamma (p371)
  gam.func <- function(x,y,mu,sig){
    u <- (x+y)/(sqrt(2))
    vee <- 1/(sig^2)
    cval <- (sqrt(2)*mu*(1-(1/vee)))^(vee-1)
    cval <- cval*exp(1-vee)
    if (vee<1) cval <- 1
    gamma <- (1/cval)*(u^(vee-1))*exp((-vee*u)/(sqrt(2)*mu))
    return(gamma)
  }
  
  #b (p371); set alpha=beta here since always want symmetrical matrix
  b.func <- function(x,y,gam){
    u <- (x+y)/(sqrt(2))
    v <- (x-y)/(sqrt(2))
    alpha<-beta<-(1-gam)/(2*gam)
    b <- (((u+v)^(alpha-1))*((u-v)^(beta-1)))/(u^(alpha+beta-2))
    return(b)
  }
  
  #smooth value
  beta.func <- function(x,y,mu,sig, gam, delta){
    gam.func(x,y,mu,sig)*b.func(x,y,gam)+delta
  }
  
  betas <- outer(X=ages.to.use,Y=ages.to.use,
                 FUN=beta.func,mu=mu,sig=sig,gam=gam,delta=delta)
  
  return(betas)
}



#Function to scale the contact matrix to a particular
#R0 given a particular population 
#
#Parameters -
#    R0 - the reproductive rate to scale to
#    state - a population at disease free equilibrium (susceptibles by age class)
#    waifw - the matrix to scale (same age class as susceptibles)
#
#Returns -
#    a scaled vertions of waifw

scale.WAIFW_base <- function(R0, DFE.state, waifw, frequency.dep=F) {
  
  if (frequency.dep) denom <- sum(DFE.state) else denom <- 1
  
  next.gen <- DFE.state*waifw/denom
  
  #get the first eigen value
  cur.R0 <- Re(eigen(next.gen)$value[1])
  
  #More correct transform
  R.ratio <- R0/cur.R0; #print(R0); #print(cur.R0); #print(R.ratio)
  waifw <- R.ratio*waifw
  
  return(list(waifw,R.ratio))
}

## Thins the polymod contact matrix for the desired age range
thin_polymod <- function(polymod, age_lower=0, age_upper=20, thin_propn=0.2,
                         same_age=TRUE){
  
  polymod1 <- polymod
  tmp_part <- polymod$participants
  tmp_contacts <- polymod$contacts
  
  tmp_contacts <- tmp_part %>% 
    select(part_id, part_age) %>% 
    full_join(tmp_contacts) %>%
    mutate(in_range=ifelse(part_age <= age_upper & part_age >= age_lower &
                             same_age & (cnt_age_exact <= age_upper & cnt_age_exact >= age_lower),
                             TRUE, FALSE))
  
  tmp_contacts_subset <- tmp_contacts %>% 
    filter(in_range == TRUE) %>%
      group_by(part_id) %>%
      sample_frac(1-thin_propn)
  
  tmp_contacts_rest <- tmp_contacts %>% filter(in_range==FALSE)
  tmp_contacts <- bind_rows(tmp_contacts_subset, tmp_contacts_rest) %>% select(-part_age)
  polymod1$contacts <- tmp_contacts
  polymod1
  
}



flatten_polymod <- function(polymod){
  
  polymod1 <- polymod
  tmp_part <- polymod$participants
  tmp_contacts <- polymod$contacts
  
  mean_contacts <- tmp_contacts %>% 
    group_by(part_id) %>% 
    summarise(contacts=n()) %>% 
    ungroup() %>% 
    summarise(n=mean(contacts)) %>% pull(n)
  
  tmp_contacts <- tmp_contacts %>% 
    group_by(part_id) %>%
    sample_n(mean_contacts, replace=TRUE)
  
  polymod1$contacts <- tmp_contacts
  polymod1
  
}


