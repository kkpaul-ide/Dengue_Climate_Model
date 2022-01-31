#' Two strain dengue transmission model with larvae compartment
#' 
#' Model function
#' 
#' @param params Named vector/dataframe of numeric values which contains 
#' the time independent input parameters for the model. 
#' @param time Integer value specifying the number of days the model is 
#' run for.
#' @param initial_conditions Named vector which contains the initial values
#' of each of the mosquito and human compartment. 
#' @param option Either difference or differential model to run.
#' @return Return a dataframe containing values of each compartment
#' at specified time points

den_model <- function(project, params, data, initial_conditions, 
                      option = "difference") {
  
  # Extract key project and data parameters
  dt <- project$dt
  time <- project$npts
  temp <- data$temp
  hum <- data$hum
  rain <- data$rain
  
  # Initialize compartments
  s_l <- rep(0,time)
  i_l1 <- rep(0,time)
  i_l2 <- rep(0,time)
  s_p <- rep(0,time)
  i_p1 <- rep(0,time)
  i_p2 <- rep(0,time)
  s_v <- rep(0,time)
  e_v1 <- rep(0,time)
  e_v2 <- rep(0,time)
  i_v1 <- rep(0,time)
  i_v2 <- rep(0,time)
  s_h <- rep(0,time)
  e_h1 <- rep(0,time)
  e_h2 <- rep(0,time)
  i_h1 <- rep(0,time)
  i_h2 <- rep(0,time)
  ci_h1 <- rep(0,time)
  ci_h2 <- rep(0,time)
  r_h1 <- rep(0,time)
  r_h2 <- rep(0,time)
  e_h12 <- rep(0,time)
  e_h21 <- rep(0,time)
  i_h12 <- rep(0,time)
  i_h21 <- rep(0,time)
  r <- rep(0,time)
  
  Nl <- rep(0,time)
  Np <- rep(0,time)
  Nv <- rep(0,time)
  Nh <- rep(0,time)
  
  newInfections <- rep(0,time)
  reInfections12 <- rep(0,time)
  reInfections21 <- rep(0,time)
  denMortality <- rep(0,time)
  
  if (option == "difference") {
    # Then run as a difference equation
    
    # Set-up parameters 
    a <- 1 - exp(-a_t(temp) * dt)
    PDR <- 1 - exp(-PDR_t(temp) * dt)
    muv <- 1 - exp(-muv_th(temp,hum,dt) * dt)
    LMR <- 1 - exp(-LMR_t(temp) * dt)
    PMR <- 1 - exp(-PMR_t(temp) * dt)
    mul_tr <- 1 - exp(-mul_tr_t(temp,rain) * dt)
    mup_tr <- 1 - exp(-mup_tr_t(temp,rain) * dt)
    EFD <- 1 - exp(-EFD_t(temp) * dt)
    
    # Birth rate and death per 1000 er year
    BR <- rep(1 - exp(-(params$BR/(1000*360))*dt),time)
    DR <- rep(1 - exp(-(params$DR/(1000*360))* dt), time)
    
    # inflow and outflow per day
    im  <- rep(params$im, time)
    em <- im
    influx <- rep(params$influx, time)
    
    # Set-up initial conditions
    s_l[1] <- initial_conditions$s_l
    i_l1[1] <- initial_conditions$i_l1
    i_l2[1] <- initial_conditions$i_l2
    s_p[1] <- initial_conditions$s_p
    i_p1[1] <- initial_conditions$i_p1
    i_p2[1] <- initial_conditions$i_p2
    s_v[1] <- initial_conditions$s_v
    e_v1[1] <- initial_conditions$e_v1
    e_v2[1] <- initial_conditions$e_v2
    i_v1[1] <- initial_conditions$i_v1
    i_v2[1] <- initial_conditions$i_v2
    s_h[1] <- initial_conditions$s_h
    e_h1[1] <- initial_conditions$e_h1
    e_h2[1] <- initial_conditions$e_h2
    i_h1[1] <- initial_conditions$i_h1
    i_h2[1] <- initial_conditions$i_h2
    ci_h1[1] <- initial_conditions$ci_h1
    ci_h2[1] <- initial_conditions$ci_h2
    r_h1[1] <- initial_conditions$r_h1
    r_h2[1] <- initial_conditions$r_h2
    e_h12[1] <- initial_conditions$e_h12
    e_h21[1] <- initial_conditions$e_h21
    i_h12[1] <- initial_conditions$i_h12
    i_h21[1] <- initial_conditions$i_h21
    r[1] <- initial_conditions$r
    
    Nl[1] <- s_l[1] + i_l1[1] + i_l2[1]
    Np[1] <- s_p[1] + i_p1[1] + i_p2[1]
    Nv[1] <- s_v[1] + e_v1[1] + e_v2[1] + i_v1[1] + i_v2[1]
    Nh[1] <- s_h[1] + e_h1[1] + e_h2[1] + i_h1[1] + i_h2[1] + ci_h1[1] +  ci_h2[1] + r_h1[1] + r_h2[1] + e_h12[1] + e_h21[1] + i_h12[1] + i_h21[1] +  r[1]
    
    for (tt in 2:time) {
      
      # Calculate time-varying non-vectorised parameters
      pMItt <- pMI(temp[tt], pMI_T0=params$pMI_T0, pMI_Tm=params$pMI_Tm, 
                   pMI_m=params$pMI_m, pMI_c=params$pMI_c)
      pIMtt <- pIM(temp[tt], pIM_T0=params$pIM_T0, pIM_Tm=params$pIM_Tm, 
                   pIM_m=params$pIM_m, pIM_c=params$pIM_c)
      FELtt <- FEL(rain[tt])
      
      # mosquito population dynamics
      s_l[tt] <- s_l[tt-1] + (EFD[tt]*FELtt*max((1-Nl[tt-1]/(Nh[tt-1]*C_t(rain[tt]))),0)*0.5*Nv[tt-1])*(1 - params$nu*(i_v1[tt-1] + i_v2[tt-1])/Nv[tt-1]) - 2*LMR[tt]*s_l[tt-1] -  mul_tr[tt]*s_l[tt-1] 
      i_l1[tt] <- i_l1[tt-1] + (EFD[tt]*FELtt*max((1-Nl[tt-1]/(Nh[tt-1]*C_t(rain[tt]))),0)*0.5*Nv[tt-1])*params$nu*i_v1[tt-1]/Nv[tt-1] - LMR[tt]*i_l1[tt-1] - mul_tr[tt]*i_l1[tt-1]
      i_l2[tt] <- i_l2[tt-1] + (EFD[tt]*FELtt*max((1-Nl[tt-1]/(Nh[tt-1]*C_t(rain[tt]))),0)*0.5*Nv[tt-1])*params$nu*i_v2[tt-1]/Nv[tt-1] - LMR[tt]*i_l2[tt-1] - mul_tr[tt]*i_l2[tt-1]
      s_p[tt] <- s_p[tt-1] + 2*LMR[tt]*s_l[tt-1] - 2*PMR[tt]*s_p[tt-1] - mup_tr[tt]*s_p[tt-1]
      i_p1[tt] <- i_p1[tt-1] + LMR[tt]*i_l1[tt-1] - PMR[tt]*i_p1[tt-1] - mup_tr[tt]*i_p1[tt-1]
      i_p2[tt] <- i_p2[tt-1] + LMR[tt]*i_l2[tt-1] - PMR[tt]*i_p2[tt-1] - mup_tr[tt]*i_p2[tt-1]
      s_v[tt] <- s_v[tt-1] + 2*PMR[tt]*s_p[tt-1] - a[tt]*pMItt*((i_h1[tt-1] + params$ADE*i_h21[tt-1])+(i_h2[tt-1] + params$ADE*i_h12[tt-1]))*s_v[tt-1]/Nh[tt-1] - muv[tt]*s_v[tt-1]
      e_v1[tt] <- e_v1[tt-1] + a[tt]*pMItt*(i_h1[tt-1] + params$ADE*i_h21[tt-1])/Nh[tt-1]*s_v[tt-1] - PDR[tt]*e_v1[tt-1] - muv[tt]*e_v1[tt-1]
      e_v2[tt] <- e_v2[tt-1] + a[tt]*pMItt*(i_h2[tt-1] + params$ADE*i_h12[tt-1])/Nh[tt-1]*s_v[tt-1] - PDR[tt]*e_v2[tt-1] - muv[tt]*e_v2[tt-1]
      i_v1[tt] <- i_v1[tt-1] + PMR[tt]*i_p1[tt-1] + PDR[tt]*e_v1[tt-1] - muv[tt]*i_v1[tt-1]
      i_v2[tt] <- i_v2[tt-1] + PMR[tt]*i_p2[tt-1] + PDR[tt]*e_v2[tt-1] - muv[tt]*i_v2[tt-1]
      
      # host population dynamics
      s_h[tt] <- s_h[tt-1] + BR[tt]*Nh[tt-1] - a[tt]*pIMtt*((i_v1[tt-1] + i_v2[tt-1])/Nh[tt-1])*s_h[tt-1] - DR[tt]*s_h[tt-1]  #+ im[tt]*Nh[tt-1] - em[tt]*s_h[tt-1] 
      e_h1[tt] <- e_h1[tt-1] + a[tt]*pIMtt*(i_v1[tt-1]/Nh[tt-1])*s_h[tt-1] - params$alpha1*e_h1[tt-1] - DR[tt]*e_h1[tt-1] + influx[tt-1] #- em[tt]*e_h1[tt-1] 
      e_h2[tt] <- e_h2[tt-1] + a[tt]*pIMtt*(i_v2[tt-1]/Nh[tt-1])*s_h[tt-1] - params$alpha2*e_h2[tt-1] - DR[tt]*e_h2[tt-1] #- em[tt]*e_h2[tt-1]
      i_h1[tt] <- i_h1[tt-1] + params$alpha1*e_h1[tt-1] - params$gamma1*i_h1[tt-1] - DR[tt]*i_h1[tt-1] #- em[tt]*i_h1[tt-1]
      i_h2[tt] <- i_h2[tt-1] + params$alpha2*e_h2[tt-1] - params$gamma2*i_h2[tt-1] - DR[tt]*i_h2[tt-1] #- em[tt]*i_h2[tt-1]
      ci_h1[tt] <- ci_h1[tt-1] + params$gamma1*i_h1[tt-1] - params$sigma1*ci_h1[tt-1] - DR[tt]*ci_h1[tt-1] #- em[tt]*ci_h1[tt-1]
      ci_h2[tt] <- ci_h2[tt-1] + params$gamma2*i_h2[tt-1] - params$sigma2*ci_h2[tt-1] - DR[tt]*ci_h2[tt-1] #- em[tt]*ci_h2[tt-1]
      r_h1[tt] <- r_h1[tt-1] + params$sigma1*ci_h1[tt-1] - a[tt]*pIMtt*(i_v2[tt-1]/Nh[tt-1])*r_h1[tt-1] - DR[tt]*r_h1[tt-1] #- em[tt]*r_h1[tt-1]
      r_h2[tt] <- r_h2[tt-1] + params$sigma2*ci_h2[tt-1] - a[tt]*pIMtt*(i_v1[tt-1]/Nh[tt-1])*r_h2[tt-1] - DR[tt]*r_h2[tt-1] #- em[tt]*r_h2[tt-1]
      e_h12[tt] <- e_h12[tt-1] + a[tt]*pIMtt*(i_v2[tt-1]/Nh[tt-1])*r_h1[tt-1] - params$alpha2*e_h12[tt-1] - DR[tt]*e_h12[tt-1] #- em[tt]*e_h12[tt-1]
      e_h21[tt] <- e_h21[tt-1] + a[tt]*pIMtt*(i_v1[tt-1]/Nh[tt-1])*r_h2[tt-1] - params$alpha1*e_h21[tt-1] - DR[tt]*e_h21[tt-1] #- em[tt]*e_h21[tt-1]
      i_h12[tt] <- i_h12[tt-1] + params$alpha2*e_h12[tt-1] - params$DIM*i_h12[tt-1] - params$gamma2*i_h12[tt-1] - DR[tt]*i_h12[tt-1] #- em[tt]*i_h12[tt-1]
      i_h21[tt] <- i_h21[tt-1] + params$alpha1*e_h21[tt-1] - params$DIM*i_h21[tt-1] - params$gamma1*i_h21[tt-1] - DR[tt]*i_h21[tt-1] #- em[tt]*i_h21[tt-1]
      r[tt] <- r[tt-1] + params$gamma2*i_h12[tt-1] + params$gamma1*i_h21[tt-1] - DR[tt]*r[tt-1] #- em[tt]*r[tt-1]
      
      # Update population size
      Nl[tt] <- s_l[tt] + i_l1[tt] + i_l2[tt]
      Np[tt] <- s_p[tt] + i_p1[tt] + i_p2[tt]
      Nv[tt] <- s_v[tt] + e_v1[tt] + e_v2[tt] + i_v1[tt] + i_v2[tt]
      Nh[tt] <- s_h[tt] + e_h1[tt] + e_h2[tt] + i_h1[tt] + i_h2[tt] + ci_h1[tt] +  ci_h2[tt] + r_h1[tt] + r_h2[tt] + e_h12[tt] + e_h21[tt] + i_h12[tt] + i_h21[tt] +  r[tt]
      
      # Extra outputs
      newInfections[tt] <- a[tt]*pIMtt*((i_v1[tt-1] + i_v2[tt-1])/Nh[tt-1])*s_h[tt-1]
      reInfections12[tt] <- a[tt]*pIMtt*(i_v2[tt-1]/Nh[tt-1])*r_h1[tt-1]
      reInfections21[tt] <- a[tt]*pIMtt*(i_v1[tt-1]/Nh[tt-1])*r_h2[tt-1]
      denMortality[tt] <-  params$DIM*i_h12[tt-1] + params$DIM*i_h21[tt-1]
      
    }
    
    # Store outputs from model in results
    return(tibble(time = 1:time, Nl, s_l, i_l1, i_l2, Np, s_p, i_p1, i_p2, Nv, s_v, e_v1, e_v2, i_v1, i_v2, 
                      Nh, s_h, e_h1, e_h2, i_h1, i_h2, ci_h1, ci_h2, 
                      r_h1, r_h2, e_h12, e_h21, i_h12, i_h21, r, newInfections,
                      reInfections12, reInfections21, denMortality))
    
  } else if (option == "ode") {
    
    source(code/den_model_ode.R)  
    
  } else {
    # Unknown option this is an error 
    stop("Unknown model equation option")
    
  }
  
  # Return the final results as a big data frame
  return(results)
  
}





