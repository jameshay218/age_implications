sigmoid_func <- function (x, k=0.1, x0=0, A=0) {
  y <- A + (1-A)/(1 + exp(-k * (x - x0)))
}


## Expand contact matrix for immunity classes
## For a given age-specific contact matrix and population size matrix (rows for age group, columns for immunity classes)
## Generates a normalised contact matrix giving contact rates between each age and immunity class combination
setup_C <- function(C1, Ns, beta_scales=rep(1, nrow=Ns)){
  Nimmunity <- ncol(Ns)
  Nage <- nrow(Ns)
  M <-  kron(C1,ones(Nimmunity,Nimmunity)) #' Non-normalised contact matrix scaled for each immunity class
  
  propns <- Ns/rowSums(Ns) #' Age/immunity population size as proportion of age population size
  propns_1 <- repmat(matrix(rep(propns, each=Nimmunity),
                            ncol = Nimmunity,byrow=FALSE),n=1,m=Nage)
  C <- M*propns_1 #' Generate scaled contact rates for age and immunity groups
  
  beta_scales_C <- repmat(matrix(rep(beta_scales, each=Nimmunity),ncol=Nimmunity*Nage,byrow=TRUE),n=Nage*Nimmunity,m=1)
  C <- C*beta_scales_C
  N_long <- c(t(Ns))
  C <- t(apply(C, 1, function(x) x/N_long))
  C[!is.finite(C)] <- 0
  #C <- t(C)
  
  return(C)
}
setup_C_explicit <- function(C1, Ns, beta_scales=rep(1, nrow=Ns)){
  Nimmunity <- ncol(Ns)
  Nage <- nrow(Ns)
  M <-  kron(C1,ones(Nimmunity,Nimmunity)) #' Non-normalised contact matrix scaled for each immunity class
  
  propns <- Ns/rowSums(Ns) #' Age/immunity population size as proportion of age population size
  propns_1 <- repmat(matrix(rep(propns, each=Nimmunity),
                            ncol = Nimmunity,byrow=FALSE),n=1,m=Nage)
  C <- M*propns_1 #' Generate scaled contact rates for age and immunity groups
  
  beta_scales_C <- repmat(matrix(rep(beta_scales, each=Nimmunity),ncol=Nimmunity*Nage,byrow=TRUE),n=Nage*Nimmunity,m=1)
  C <- C*beta_scales_C
  return(C)
}
## Get's the transmissibility parameter required for a desired R0
get_beta <- function(C, propns, inf_period=5, r0){
  fi <- kron(propns, ones(1, length(propns)))
  fj <- kron(t(propns), ones(length(propns),1))
  M <- C * fi/fj
  R0 <- max(eigen(M)$values)
  beta <- r0/(R0*inf_period)
  beta
}

## Get's the transmissibility parameter required for a desired R0
get_beta_vector <- function(C, propns, inf_period=5, r0, beta_scales){
  fi <- kron(propns, ones(1, length(propns)))
  fj <- kron(t(propns), ones(length(propns),1))
  C1 <- t(beta_scales * t(C))
  M <- C1 * fi/fj
  R0 <- max(eigen(M)$values)
  beta <- r0/(R0*inf_period)
  beta
}

get_sigmoid_func_outputs <- function(ages, pars){
  final <- matrix(nrow=nrow(pars),ncol=length(ages))
  for(j in 1:nrow(pars)){
    A <- pars$minimum_transmissibility_scale[j]
    k <- pars$transmissibility_scale[j]
    x0 <- pars$transmissibility_midpoint[j]
    y <- sigmoid_func(ages, k, x0, A)
    final[j,] <- y
  }
  colnames(final) <- ages
  final <- cbind(pars, final)
  final <- final %>% pivot_longer(as.character(ages),names_to="age_group",values_to="beta")
  return(final)
}

## Thins the polymod contact matrix for the desired age range
thin_polymod <- function(polymod, age_lower=0, age_upper=20, thin_propn=0.2,
                         same_age=TRUE){
  
  polymod1 <- polymod
  tmp_part <- polymod$participants
  tmp_contacts <- polymod$contacts
  tmp_contacts <- tmp_part %>% 
    select(part_id, part_age) %>% 
    right_join(tmp_contacts) %>% 
    mutate(row_id=1:n())
  ## Get contacts with age within range or overlapping range
  
  if(same_age){
    tmp_contacts_tothin <- tmp_contacts %>%
      filter((part_age <= age_upper & part_age >= age_lower) & 
               same_age & ((cnt_age_exact <= age_upper & cnt_age_exact >= age_lower) |
                           (cnt_age_est_min < age_upper & cnt_age_est_min > age_lower) |
                           (cnt_age_est_max > age_lower & cnt_age_est_max < age_upper)))
  } else {
    tmp_contacts_tothin <- tmp_contacts %>%
      filter((part_age <= age_upper & part_age >= age_lower))
  }
  tmp_contacts_subset <- tmp_contacts_tothin %>% 
    group_by(part_id) %>%
    sample_frac(1-thin_propn)
  
  tmp_contacts_rest <- tmp_contacts %>% filter(!(row_id %in% tmp_contacts_tothin$row_id))
  tmp_contacts <- bind_rows(tmp_contacts_subset, tmp_contacts_rest) %>% select(-part_age, -row_id)
  polymod1$contacts <- tmp_contacts
  polymod1
}

