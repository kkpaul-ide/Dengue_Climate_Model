#' General function for the Briere fit
#' 
#' @param x temperature at which mosquito trait will be estimated
#' @param c rate constant
#' @param T0 critical thermal minimum
#' @param Tm critical thermal maximum
briere <-  function(t, tmin, tmax, m, c) {
        # Use ifelse to vectorise
        return(ifelse(t < tmin | t > tmax, 0, c*t*(t-tmin)*(tmax-t)^(1/m)))
}


#' General function for the inverted quadratic fit.
#' 
#' @param x temperature at which mosquito trait will be estimated
#' @param c rate constant
#' @param T0 critical thermal minimum
#' @param Tm critical thermal maximum
inverted_quadratic <- function(x, c, T0, Tm, timestep){
        return(ifelse(x < T0 | x > Tm, 1.0/timestep, 1.0/(c*(x-T0)*(x-Tm))))
}


#' This is the general function for the quadratic fit.
#'
#' @param x temperature at which mosquito trait will be estimated
#' @param c rate constant
#' @param T0 critical thermal minimum
#' @param Tm critical thermal maximum 
quadratic <- function(x, c, T0, Tm){
        if((x < T0) | (x > Tm))
                0.0
        else
                c*(x-T0)*(x-Tm)
}

#' General function for the polynomial fit
#' 
#' @param x temperature at which mosquito trait will be estimated
#' @param b0 constant term
#' @param b1 first coefficient
#' @param b2 second coefficient
#' @param b3 third coefficient
#' @param b4 fourth coefficient
#' @param b5 fifth coefficient
#' @param b6 sixth coefficient
#' @param b7 seventh coefficient
#' @param b8 eighth coefficient
poly <- function(x, b0=0, b1=0, b2=0, b3=0, b4=0, b5=0, b6=0, b7=0, b8=0){
        return(b0 + b1*x + b2*x^2 + b3*x^3 + b4*x^4 + b5*x^5 + b6*x^6 + b7*x^7 + b8*x^8)  
}

#' function to convert rate to probability
prob <- function(x, timestep){
        1-exp(-x * timestep)
}

#' Convert parameters into a list - Might be a bit convoluted but I prefer 
#' things in columns 
params2df <- function(pm, columnName) {
        return(pm %>% select(parameter, columnName) %>% spread(1,2))
}


#' Fill values for selected time points in a day
#' 
#' @param x vector containing daily climate data
#' @param ts timestep in a day for which values need to created
#' @return vector with filled gaps for each timestep
# fillGap <- function(x,ts){
#         position <- c(2:(length(x)))
#         anyList <- list()
#         for(i in 1:(length(x)-1)){
#                 if (i < length(x)){
#                         anyList[[i]] <- seq(x[i], x[i+1], length.out = 1/ts-1)
#                 } else {
#                         anyList[[i]] <- seq(x[i], x[i], length.out = 1/ts-1)
#                 }
#         }
#         insert(x, at = position, values = anyList)
# }

fillGap <- function(clim_var,t_step,dtr,option){
                if (option=="temp"){
                sinusoidal_temp <- function (x,y,z){
                        return(x + y*sin(0:((1/z)-1)*pi*z*2)/2)
                }
                aa <- list()
                for (ii in 1:length(clim_var)){
                        if (ii < length(clim_var)) {
                        aa[[ii]] <- sinusoidal_temp(clim_var[ii],dtr[ii],t_step)
                        } else {
                                aa[[ii]] <- clim_var[ii]
                        }
                }
                return(unlist(aa))
                
        } else if (option=="hum"){
                
                return(c(rep(clim_var[1:length(clim_var)-1],each=1/t_step),clim_var[length(clim_var)]))
        }  else if (option=="rain") {
                return(c(rep((clim_var*t_step)[1:length(clim_var)-1],each=1/t_step),(clim_var*t_step)[length(clim_var)]))
        }
}



#' Mosquito infection per bite on an infectious host ----------------------
#' 
#' briere fit to empirical data #### REFERENCE! WHICH MODEL!
#' 
#' 
#' 
#' function of temperature
#' @param t daily mean temperature
#' @return return probability of mosquito infection per bite of an infectious host
pMI <- function(t, pMI_T0=16, pMI_Tm=36.6, pMI_m=3.3, pMI_c=9.8e-04){
        return(briere(t, pMI_T0, pMI_Tm, pMI_m, pMI_c))
}

#' Human infection per bite on an infectious mosquito
#' 
#' briere fit to empirical data
#' 
#' function of temperature
#' @param t daily mean temperature
#' @return return probability of human infection per bite of an infectious mosquito

pIM <- function(t, pIM_T0=18.6, pIM_Tm=39, pIM_m=2.03, pIM_c=7.03e-04){
        return(briere(t, pIM_T0, pIM_Tm, pIM_m, pIM_c))
}


#' Daily biting rate
#' 
#' briere fit to empirical data
#' 
#' function of temperature
#' @param t daily mean temperature
#' @return return mosquito biting rate per day
a_t <- function(t, a_T0=15.9, a_Tm=40.1, a_m=1.11, a_c=8.7e-5){
        return(briere(t, a_T0, a_Tm, a_m, a_c))
}

#' Larvae maturation rate
#' 
#' Briere fit to empirical data
#' 
#' function of temperature
#' @param t daily mean temperature
#' @return return rate of larvae turn into pupae per day
LMR_t <- function(t, LMR_T0=10, LMR_Tm=40, LMR_m=2, LMR_c=11.96e-5) {
        return(briere(t, LMR_T0, LMR_Tm, LMR_m, LMR_c))
}


#' Pupae maturation rate
#' 
#' Briere fit to empirical data
#' 
#' function of temperature
#' @param t daily mean temperature
#' @return return rate of pupae turn into adult per day
PMR_t <- function(t, PMR_T0=10, PMR_Tm=40, PMR_m=3.31, PMR_c=5e-4) {
        return(briere(t, PMR_T0, PMR_Tm, PMR_m, PMR_c))
}


#' Adult mosquito mortality rate
#' 
#' fitted spline model based on a pooled survival analysis of Ae. aegypti
#' 
#' function of temperature and humidity
#' @param temp daily mean temperature
#' @param hum daily saturation vapour pressure deficit calculated from mean temperature
#' and relative humidity
#' @return return daily adult mosquito mortality rate
#' from Caldwell et al

muv_th <- function(temp, hum, timestep){
        return(ifelse(hum <= 1,inverted_quadratic(temp,-1.48e-01, 9.16, 37.73, timestep)+(1-(0.01256 + 2.00893*hum))*0.005,
                      inverted_quadratic(temp,-1.48e-01, 9.16, 37.73, timestep)+(1-(1.2248 + 0.2673*hum))*0.01))
}

#' virus incubation rate/Parasite development rate, inverse of extrinsic incubation period
#' 
#' briere fit to empirical data
#' 
#' function of temperature
#' @param t daily mean temperature
#' @return return daily rate of mosquito being infectious
#' Mordecai et al. (2017)

d_VIR <- function(t,  VIR_T0=16, VIR_Tm=45, VIR_m =1.91, VIR_c=8.6e-05){
        return(briere(t, VIR_T0, VIR_Tm, VIR_m, VIR_c))
}

#' parasite development rate, inverse of extrinsic incubation period
#' 
#' Briere fit to empirical data
#' 
#' function of temperature
#' @param t daily mean temperature
#' @return return daily rate of mosquito being infectious
#' 
#' Tjaden et al. (2013)
# PDR_t <- function(t,PDR_b0=-1.678, PDR_b1=0.344, PDR_b2=-2.422e-2,PDR_b3=7.252e-4,PDR_b4=-7.713e-6){
#   
#   return(poly(t, PDR_b0, PDR_b1, PDR_b2, PDR_b3, PDR_b4))
#   
# }

PDR_t <- function(t, PDR_T0=10.68, PDR_Tm=45.90, PDR_m = 2, PDR_c=6.56e-05){
        return(briere(t, PDR_T0, PDR_Tm, PDR_m, PDR_c))
}

#' Larvae mortality rate
#' 
#' From Lee et al ....
#' 
#' function of rain
#' @param pr amount of rainfall summed over the prior two-week period
#' @param n total number of larvae
#' Note: Typos in Lee et al.(2018) 
#' pr_max and pr_min must be assigned

d_mul_tr <- function(t,pr,n,pr_min,pr_max) {
        mul <- 0.08 # minimum larvae mortality rate
        K0 <- 250000 # standard carrying capacity
        P_norm <- (pr-pr_min)/(pr_max-pr_min)
        return((exp(-t/2)+mul)*(1+n/K0*(P_norm +1)))
}

#' Larvae mortality rate
#' 
#' function of temperature and rain
#' @param t daily mean temperature
#' @param pr amount of rainfall summed over the prior two-week period
#' Note: anonymous function theta and temperature dependent mortality included
mul_tr_t <- function(t,pr,pr_c=30,mu_al=0.001,mul_b0=2.32,mul_b1=-4.19e-1,mul_b2=2.73e-2,mul_b3=-7.53e-4,mul_b4=7.50e-6){
        theta <- function(r,r_c=30){
                return(ifelse(r - r_c > 0 | r - r_c == 0, 1, 0))
        }
        mul_t <- function(x){
                return(poly(t, mul_b0, mul_b1, mul_b2, mul_b3, mul_b4))
        }
        return(mul_t(t)*(1+mu_al*(pr-pr_c)*theta(pr)))
}

#' Pupae mortality rate
#' 
#' Polynomial fit to empirical data
#' 
#' function of temperature
#' @param t daily mean temperature
#' @return return daily pupae mortality rate
#' Note: Liu Helmersson (2019) shows mup_b3 to be 4.39e-7
mup_tr_t <- function(t,pr,pr_c=30,mu_ap=0.001,mup_b0=4.25e-1,mup_b1=-3.25e-2,mup_b2=7.06e-4,mup_b3=4.39e-7){
        theta <- function(r,r_c=30){
                return(ifelse(r - r_c > 0 | r - r_c == 0, 1, 0))
        }
        mup_t <- function(x) {
                return(poly(t, mup_b0, mup_b1, mup_b2, mup_b3))
        }
        return(mup_t(t)*(1+mu_ap*(pr-pr_c)*theta(pr)))
}

#' Oviposition rate/eggs laid per female mosquito
#' 
#' polynomial fit to empirical data
#' @param t daily mean temperature
#' @return return Oviposition rate per female mosquito
# EFD_t <- function(t,EFD_b0=-5.40, EFD_b1=1.80, EFD_b2=-2.12e-1,EFD_b3=1.02e-2,EFD_b4=-1.51e-4){
#         return(poly(t, EFD_b0, EFD_b1, EFD_b2, EFD_b3, EFD_b4))
# }

#'Briere fit to empirical data
#' 
#' function of temperature
#' @param t daily mean temperature
#' @return return daily rate of eggs laid per female
EFD_t <- function(t, EFD_T0=16.1, EFD_Tm=34.4, EFD_m=1.76, EFD_c=9e-3) {
        return(briere(t, EFD_T0, EFD_Tm, EFD_m, EFD_c))
}


#' Fraction of eggs hatching to larvae, function of rain
#' 
#' 
#' @param pr daily total rainfall
#' @return fraction of eggs hatching to larvae
FEL <- function(pr,q0=0.2,q1=0.02,q2=0.037){
        return((q1*pr/(q0+q1*pr))+q2)
        
}

#' larval carrying capacity per breeding site
#' 
#' typo in Liu Helmersson et al., Correct form available in Yang et al. (2016)
#' 
#' 
C_t <- function(pr,c0=5,c1=30,c2=0.1){
        return((c0*pr/(c1+pr))+c2)
}