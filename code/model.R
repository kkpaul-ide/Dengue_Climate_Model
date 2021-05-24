#' Two strain dengue transmission model
#' 
#' Model function
#' 
#' @param params Named vector/dataframe of numeric values which contains 
#' the time independent input parameters for the model. 
#' @param time Integer value specifying the number of days the model is 
#' run for.
#' @param intial_conditions Named vector which contains the initial values
#' of each of the mosquito and human compartment. 
#' @param option Either difference or differential model to run.
#' @return Return a dataframe containing values of each compartment
#' at specified time points



model <- function(params, time, intial_conditions, option = "difference") {
  
  # Initalize results
  # results <- data.frame(s = rep(0,length(time)))
  
  if (option == "difference") {
    # Then run as a difference equation with daily time step
    
    s_l[1] <- intial_conditions$s_l
    i_l1[1] <- intial_conditions$i_l1
    i_l2[1] <- intial_conditions$i_l2
    s_v[1] <- intial_conditions$s_v
    e_v1[1] <- intial_conditions$e_v1
    e_v2[1] <- intial_conditions$e_v2
    i_v1[1] <- intial_conditions$i_v1
    i_v2[1] <- intial_conditions$i_v2
    s_h[1] <- intial_conditions$s_h
    e_h1[1] <- intial_conditions$e_h1
    e_h2[1] <- intial_conditions$e_h2
    i_h1[1] <- intial_conditions$i_h1
    i_h2[1] <- intial_conditions$i_h2
    ci_h1[1] <- intial_conditions$ci_h1
    ci_h2[1] <- intial_conditions$ci_h2
    r_h1[1] <- intial_conditions$r_h1
    r_h2[1] <- intial_conditions$r_h2
    e_h12[1] <- intial_conditions$e_h12
    e_h21[1] <- intial_conditions$e_h21
    i_h12[1] <- intial_conditions$i_h12
    i_h21[1] <- intial_conditions$i_h21
    r[1] <- intial_conditions$r
    
    
    for (t in 1:length(time)) {
      
      # General function for the Briere fit.
      briere <- function(x, c, T0, Tm){
        if((x < T0) | (x > Tm))
          0.0
        else
          c*x*(x-T0)*sqrt(Tm-x)
      }
      
      # General function for the inverted quadratic fit.
      inverted_quadratic <- function(x, c, T0, Tm, timestep){
        if((x < T0) | (x > Tm))
          1.0/timestep
        else
          1.0/(c*(x-T0)*(x-Tm))
      }
      
      # Model parameters affected by temperature, humidity, and rainfall over time
      # Adult mosquito mortality rate (1/adult lifespan)
      # depend on temperature and humidity
      muv_th <- function(temp, hum, timestep){
        if (hum <= 1){
          inverted_quadratic(temp,-1.48e-01, 9.16, 37.73, timestep)+
            (1-(0.01256 + 2.00893*hum))*0.005
        } else {
          inverted_quadratic(temp,-1.48e-01, 9.16, 37.73, timestep)+
            (1-(1.2248 + 0.2673*hum))*0.01
        }
      }
      # Rate larvae turn into adults per day, pre-adult maturation rate, function
      # of temperature
      PMR <- function(temp){
        briere(temp,7.86e-05,11.36,39.17)
      }
      # biting rate
      a <- function(temp){
        briere(temp,2.02e-04,13.35,40.08)
      }
      # probability of mosquito infection per	bite on	an infectious host
      pMI <- function(temp){
        briere(temp,4.91e-04,12.22,37.46)
      }
      # probability of human infection per bite by an	infectious mosquito
      pIM <- function(temp){
        briere(temp,8.49e-04,17.05,35.83)
      }
      # virus incubation rate
      VIR <- function(temp){
        briere(temp,6.56e-05,10.68,45.90)
      }
      # Lambda calculations - these are time dependent and incorporate the 
      # effects of temperature on vector to host transmission. 
      # Beta equations as well
      lambda1 <- a(temp[t])*pIM(temp[t])*iv_1[t]/Nh[t]
      lambda2 <- a(temp[t])*pIM(temp[t])*iv_1[t]/Nh[t]
      
      # host population dynamics
      s_h[t]   <- s_h[t-1]   + params$mu_hb*Nh[t] - (lambda1 + lambda2)*s_h[t-1] -
        params$mu_h*s_h[t-1]  
      e_h1[t]  <- e_h1[t-1]  + lambda1*s_h[t-1]   - (alpha1 + mu_h)*e_h1[t-1]
      e_h2[t]  <- e_h2[t-1]  + lambda2*s_h[t-1]   - (alpha2 + mu_h)*e_h2[t-1]
      i_h1[t]  <- i_h1[t-1]  + alpha1*e_h1[t-1]   - (gamma1 + mu_h)*i_h1[t-1]
      i_h1[t]  <- i_h2[t-1]  + alpha2*e_h2[t-1]   - (gamma2 + mu_h)*i_h2[t-1]
      ci_h1[t] <- ci_h1[t-1] + gamma1*i_h1[t-1]   - (sigma1 + mu_h)*ci_h1[t-1]
      ci_h2[t] <- ci_h2[t-1] + gamma2*i_h2[t-1]   - (sigma2 + mu_h)*ci_h2[t-1]
      r_h1[t]  <- r_h1[t-1]  + sigma1*ci_h1[t-1]  - lambda2*r_h1[t-1] - 
        mu_h*r_h1[t-1]
      r_h2[t]  <- r_h2[t-1]  + sigma2*ci_h2[t-1]  - lambda1*r_h2[t-1] - 
        mu_h*r_h2[t-1]
      e_h12[t] <- 
      e_h21[t] <- 
      i_h12[t] <-
      i_h21[t] <-
      r[t]     <-   
      
      # mosquito population dynamics
      s_l[t]   <- s_l[t-1] 
      i_l1[t]  <- i_l1[t-1]
      i_l2[t]  <- i_l2[t-1]
      s_v[t]   <- s_v[t-1]
      e_v1[t]  <- e_v1[t-1]
      e_v2[t]  <- e_v2[t-1]
      i_v1[t]  <- i_v1[t-1]
      i_v2[t]  <- i_v2[t-1]  
        
        
        # Extra outputs
        newInfections[t] <- lambda*s[t]
      
      
    }
    
    # Store outputs from model in results
    results$s <- s
    results$new_infections <- new_infections
    
  } else if (option == "ode") {
    
    
    ds <- s 
    
    
  } else {
    # Unknown option this is an error 
    stop("Unknown model equation option")
    
  }
  
  
  
  # Return the final results as a big data frame
  return(results)
  
}
