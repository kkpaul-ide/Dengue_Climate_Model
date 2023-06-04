#' Two strain dengue transmission model with larvae compartment
#' 
#' The function extracts key project and data parameters from the input, and 
#' initializes a number of compartments to represent different stages of infection 
#' and recovery in the population. It then calculates time-varying parameters based on 
#' temperature, humidity, and rain data and model parameters, and uses these 
#' parameters to simulate the spread of the disease through the population over time.
#' 
#' @param project a list containing the parameters for the simulation, 
#' including the time step and the number of time points  
#' @param params a list of model parameters, including temperature-dependent parameters 
#' for different rates of infection, mortality, and recovery 
#' @param data a list containing temperature and humidity data, as well as data on rain
#' @param initial_conditions a list of the initial conditions for the compartments in the model
#' of each of the mosquito and human compartment. 
#' @param option a string indicating whether the model should be run as a difference 
#' equation (default) or a differential equation
#' @return a number of vectors representing the number of individuals in each compartment at each time 
#' point, as well as other summary statistics such as new infections, re-infections, and mortality
#' 
#' @author Kishor K. Paul \email{Kpaul@kirby.unsw.edu.au}
#' @author Richard T. Gray \email{Rgray@kirby.unsw.edu.au}
#' 

den_model <- function(project, params, data, initial_conditions, option = "difference") {
        
        # Extract key project and data parameters
        dt <- project$dt
        time <- project$npts
        temp <- data$temp
        hum <- data$hum
        # rain <- data$rain
        rain <-  data$rain14/14 #rep(0, 175297)#
        
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
        newInfections1 <- rep(0,time)
        newInfections2 <- rep(0,time)
        reInfections12 <- rep(0,time)
        reInfections21 <- rep(0,time)
        denMortality <- rep(0,time)
        
        if (option == "difference") {
                # Then run as a difference equation
                
                # Calculate time-varying vectorized parameters
                # a, PDR, LMR, PMR, EFD require conversion to daily probability
                a <- 1 - exp(-briere(temp, params$a_T0, params$a_Tm, params$a_m, params$a_c) * dt)
                PDR <- 1 - exp(-(briere(temp, params$PDR_T0, params$PDR_Tm, params$PDR_m, params$PDR_c)) * dt)
                LMR <- 1 - exp(-briere(temp, params$LMR_T0, params$LMR_Tm, params$LMR_m, params$LMR_c) * dt)
                PMR <- 1 - exp(-briere(temp, params$PMR_T0, params$PMR_Tm, params$PMR_m, params$PMR_c) * dt)
                EFD <- 1 - exp(-briere(temp, params$EFD_T0, params$EFD_Tm, params$EFD_m, params$EFD_c) * dt)
                
                pMI <- params$betaMI * briere(temp, params$pMI_T0, params$pMI_Tm, params$pMI_m, params$pMI_c)
                pIM <- params$betaIM * briere(temp, params$pIM_T0, params$pIM_Tm, params$pIM_m, params$pIM_c)
                
                # muv <- 1 - exp(-muv_th(temp,hum) * dt)
                # mul_tr <- 1 - exp(-mul_tr_t(temp,rain) * dt)
                # mup_tr <- 1 - exp(-mup_tr_t(temp,rain) * dt)
                #when using quadratic fit
                muv <- 1 - exp(-muv_th(temp, params$muv_c, params$muv_T0, params$muv_Tm, hum) * dt)
                mul_tr <- 1 - exp(-mul_tr_t(temp, params$mul_tr_c, params$mul_tr_T0, params$mul_tr_Tm, rain) * dt)
                mup_tr <- 1 - exp(-mup_tr_t(temp, params$mup_tr_c, params$mup_tr_T0, params$mup_tr_Tm, rain) * dt)
                
                FEL <- FEL(rain)
                
                # Birth and death rate per 1000 per timestep
                
                BR <- rep(1 - exp(-(project$birthRate/(1000*360))* project$dt),each=(project$npts/length(project$birthRate)))
                DR <- rep(1 - exp(-(project$deathRate/(1000*360))* project$dt),each=(project$npts/length(project$deathRate)))
                
                # inflow and outflow per day - adjust influx for start 
                # time of first and second strain. 
                im  <- rep(params$im, time)
                em <- rep(params$em, time)
                influx1 <- rep(params$influx1 * dt, time)
                influx2 <- rep(params$influx2 * dt, time)
                
                influx1[project$pts_times < lubridate::as_datetime(paste0(params$start1, "-01-01"))] <- 0
                influx2[project$pts_times < lubridate::as_datetime(paste0(params$start2, "-01-01"))] <- 0
                influx1[project$pts_times == lubridate::as_datetime(paste0(params$start1, "-01-01"))] <- 1
                influx2[project$pts_times == lubridate::as_datetime(paste0(params$start2, "-01-01"))] <- 1
                
                # Carrying capacity multiplication factor
                kf <- params$kfactor
                
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
                Nh[1] <- s_h[1] + e_h1[1] + e_h2[1] + i_h1[1] + i_h2[1] + ci_h1[1] + ci_h2[1] + r_h1[1] + r_h2[1] + e_h12[1] + e_h21[1] + i_h12[1] + i_h21[1] +  r[1]
                
                for (tt in 2:time) {
                        
                        # mosquito population dynamics
                        s_l[tt] <- s_l[tt-1] + (EFD[tt]*FEL[tt]*max((1-Nl[tt-1]/(Nh[tt-1]*kf*K_t(rain[tt]))),0)*0.5*Nv[tt-1])*(1 - params$nu*(i_v1[tt-1] + i_v2[tt-1])/ifelse(Nv[tt-1]==0,1,Nv[tt-1])) - 2*LMR[tt]*s_l[tt-1] -  mul_tr[tt]*s_l[tt-1] 
                        i_l1[tt] <- i_l1[tt-1] + (EFD[tt]*FEL[tt]*max((1-Nl[tt-1]/(Nh[tt-1]*kf*K_t(rain[tt]))),0)*0.5*Nv[tt-1])*params$nu*i_v1[tt-1]/ifelse(Nv[tt-1]==0,1,Nv[tt-1]) - LMR[tt]*i_l1[tt-1] - mul_tr[tt]*i_l1[tt-1]
                        i_l2[tt] <- i_l2[tt-1] + (EFD[tt]*FEL[tt]*max((1-Nl[tt-1]/(Nh[tt-1]*kf*K_t(rain[tt]))),0)*0.5*Nv[tt-1])*params$nu*i_v2[tt-1]/ifelse(Nv[tt-1]==0,1,Nv[tt-1]) - LMR[tt]*i_l2[tt-1] - mul_tr[tt]*i_l2[tt-1]
                        s_p[tt] <- s_p[tt-1] + LMR[tt]*s_l[tt-1] - PMR[tt]*s_p[tt-1] - mup_tr[tt]*s_p[tt-1]
                        i_p1[tt] <- i_p1[tt-1] + LMR[tt]*i_l1[tt-1] - PMR[tt]*i_p1[tt-1] - mup_tr[tt]*i_p1[tt-1]
                        i_p2[tt] <- i_p2[tt-1] + LMR[tt]*i_l2[tt-1] - PMR[tt]*i_p2[tt-1] - mup_tr[tt]*i_p2[tt-1]
                        s_v[tt] <- s_v[tt-1] + PMR[tt]*s_p[tt-1] - a[tt]*pMI[tt]*((i_h1[tt-1] + params$ADE*i_h21[tt-1])+(i_h2[tt-1] + params$ADE*i_h12[tt-1]))*s_v[tt-1]/Nh[tt-1] - muv[tt]*s_v[tt-1]
                        e_v1[tt] <- e_v1[tt-1] + a[tt]*pMI[tt]*(i_h1[tt-1] + params$ADE*i_h21[tt-1])/Nh[tt-1]*s_v[tt-1] - PDR[tt]*e_v1[tt-1] - muv[tt]*e_v1[tt-1]
                        e_v2[tt] <- e_v2[tt-1] + a[tt]*pMI[tt]*(i_h2[tt-1] + params$ADE*i_h12[tt-1])/Nh[tt-1]*s_v[tt-1] - PDR[tt]*e_v2[tt-1] - muv[tt]*e_v2[tt-1]
                        i_v1[tt] <- i_v1[tt-1] + PMR[tt]*i_p1[tt-1] + PDR[tt]*e_v1[tt-1] - muv[tt]*i_v1[tt-1]
                        i_v2[tt] <- i_v2[tt-1] + PMR[tt]*i_p2[tt-1] + PDR[tt]*e_v2[tt-1] - muv[tt]*i_v2[tt-1]
                        
                        # host population dynamics
                        s_h[tt] <- s_h[tt-1] +  BR[tt]*Nh[tt-1] - # newpop[tt] -
                                a[tt]*pIM[tt]*((i_v1[tt-1] + i_v2[tt-1])/Nh[tt-1])*s_h[tt-1] - DR[tt]*s_h[tt-1]  + im[tt]*Nh[tt-1]*(1-params$exseroprev*params$serofactor) - em[tt]*s_h[tt-1] 
                        e_h1[tt] <- e_h1[tt-1] + a[tt]*pIM[tt]*(i_v1[tt-1]/Nh[tt-1])*s_h[tt-1] - (params$alpha1 + DR[tt] + em[tt])*e_h1[tt-1]  + influx1[tt-1]  
                        e_h2[tt] <- e_h2[tt-1] + a[tt]*pIM[tt]*(i_v2[tt-1]/Nh[tt-1])*s_h[tt-1] - params$alpha2*e_h2[tt-1] - DR[tt]*e_h2[tt-1] + influx2[tt-1] - em[tt]*e_h2[tt-1]
                        i_h1[tt] <- i_h1[tt-1] + params$alpha1*e_h1[tt-1] - params$gamma1*i_h1[tt-1] - DR[tt]*i_h1[tt-1] - em[tt]*i_h1[tt-1]
                        i_h2[tt] <- i_h2[tt-1] + params$alpha2*e_h2[tt-1] - params$gamma2*i_h2[tt-1] - DR[tt]*i_h2[tt-1] - em[tt]*i_h2[tt-1]
                        ci_h1[tt] <- ci_h1[tt-1] + params$gamma1*i_h1[tt-1] - params$sigma1*ci_h1[tt-1] - DR[tt]*ci_h1[tt-1] - em[tt]*ci_h1[tt-1]
                        ci_h2[tt] <- ci_h2[tt-1] + params$gamma2*i_h2[tt-1] - params$sigma2*ci_h2[tt-1] - DR[tt]*ci_h2[tt-1] - em[tt]*ci_h2[tt-1]
                        r_h1[tt] <- r_h1[tt-1] + params$sigma1*ci_h1[tt-1] - a[tt]*pIM[tt]*(i_v2[tt-1]/Nh[tt-1])*r_h1[tt-1] - DR[tt]*r_h1[tt-1] + im[tt]*Nh[tt-1]*params$exseroprev*params$serofactor*(r_h1[tt-1]/(r_h1[tt-1]+r_h2[tt-1]+r[tt-1])) - em[tt]*r_h1[tt-1]
                        r_h2[tt] <- r_h2[tt-1] + params$sigma2*ci_h2[tt-1] - a[tt]*pIM[tt]*(i_v1[tt-1]/Nh[tt-1])*r_h2[tt-1] - DR[tt]*r_h2[tt-1] + im[tt]*Nh[tt-1]*params$exseroprev*params$serofactor*(r_h2[tt-1]/(r_h1[tt-1]+r_h2[tt-1]+r[tt-1])) - em[tt]*r_h2[tt-1]
                        e_h12[tt] <- e_h12[tt-1] + a[tt]*pIM[tt]*(i_v2[tt-1]/Nh[tt-1])*r_h1[tt-1] - params$alpha2*e_h12[tt-1] - DR[tt]*e_h12[tt-1] - em[tt]*e_h12[tt-1]
                        e_h21[tt] <- e_h21[tt-1] + a[tt]*pIM[tt]*(i_v1[tt-1]/Nh[tt-1])*r_h2[tt-1] - params$alpha1*e_h21[tt-1] - DR[tt]*e_h21[tt-1] - em[tt]*e_h21[tt-1]
                        i_h12[tt] <- i_h12[tt-1] + params$alpha2*e_h12[tt-1] - params$DIM*i_h12[tt-1] - params$gamma2*i_h12[tt-1] - DR[tt]*i_h12[tt-1] - em[tt]*i_h12[tt-1]
                        i_h21[tt] <- i_h21[tt-1] + params$alpha1*e_h21[tt-1] - params$DIM*i_h21[tt-1] - params$gamma1*i_h21[tt-1] - DR[tt]*i_h21[tt-1] - em[tt]*i_h21[tt-1]
                        r[tt] <- r[tt-1] + params$gamma2*i_h12[tt-1] + params$gamma1*i_h21[tt-1] - DR[tt]*r[tt-1] + im[tt]*Nh[tt-1]*params$exseroprev*params$serofactor*(r[tt-1]/(r_h1[tt-1]+r_h2[tt-1]+r[tt-1])) - em[tt]*r[tt-1]
                        
                        # Update population size
                        Nl[tt] <- s_l[tt] + i_l1[tt] + i_l2[tt]
                        Np[tt] <- s_p[tt] + i_p1[tt] + i_p2[tt]
                        Nv[tt] <- s_v[tt] + e_v1[tt] + e_v2[tt] + i_v1[tt] + i_v2[tt]
                        Nh[tt] <- s_h[tt] + e_h1[tt] + e_h2[tt] + i_h1[tt] + i_h2[tt] + ci_h1[tt] +  ci_h2[tt] + r_h1[tt] + r_h2[tt] + e_h12[tt] + e_h21[tt] + i_h12[tt] + i_h21[tt] +  r[tt]
                        
                        # Extra outputs
                        newInfections[tt] <- a[tt]*pIM[tt]*((i_v1[tt-1] + i_v2[tt-1])/Nh[tt-1])*s_h[tt-1]
                        newInfections1[tt] <- a[tt]*pIM[tt]*(i_v1[tt-1]/Nh[tt-1])*s_h[tt-1]
                        newInfections2[tt] <- a[tt]*pIM[tt]*(i_v2[tt-1]/Nh[tt-1])*s_h[tt-1]
                        reInfections12[tt] <- a[tt]*pIM[tt]*(i_v2[tt-1]/Nh[tt-1])*r_h1[tt-1]
                        reInfections21[tt] <- a[tt]*pIM[tt]*(i_v1[tt-1]/Nh[tt-1])*r_h2[tt-1]
                        denMortality[tt] <-  params$DIM*i_h12[tt-1] + params$DIM*i_h21[tt-1]
                        
                        # allDeaths[tt] <- 0
                        # newPop[tt] <- 0
                        
                }
                
                # Store outputs from model in results
                return(tibble(time = 1:time, Nl, s_l, i_l1, i_l2, Np, s_p, i_p1, i_p2, Nv, s_v, e_v1, e_v2, i_v1, i_v2, 
                              Nh, s_h, e_h1, e_h2, i_h1, i_h2, ci_h1, ci_h2, 
                              r_h1, r_h2, e_h12, e_h21, i_h12, i_h21, r, newInfections, newInfections1, newInfections2,
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





