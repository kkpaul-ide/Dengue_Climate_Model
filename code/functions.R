#' General function for the Briere fit
#' 
#' @param x temperature at which mosquito trait will be estimated
#' @param c rate constant
#' @param T0 critical thermal minimum
#' @param Tm critical thermal maximum
briere <-  function(t, tmin, tmax, m, c) {
        # Use ifelse to vectorise
        return(ifelse(t < tmin | t > tmax, 0, c * t * (t-tmin) * (tmax-t)^(1/m)))
}


#' General function for the inverted quadratic fit.
#' 
#' @param x temperature at which mosquito trait will be estimated
#' @param c rate constant
#' @param T0 critical thermal minimum
#' @param Tm critical thermal maximum
# inverted_quadratic <- function(x, c, T0, Tm, timestep){
#         return(ifelse(x < T0 | x > Tm, 1.0/timestep, 1.0/(c*(x-T0)*(x-Tm))))
# }

inverted_quadratic <- function(x, c, T0, Tm){
        return(ifelse(x < T0 | x > Tm, 1.0, 1.0/(c*(x-T0)*(x-Tm))))
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

# muv_th <- function(temp, hum){
#         return(ifelse(hum <= 1,inverted_quadratic(temp,-1.48e-01, 9.16, 37.73)+(1-(0.01256 + 2.00893*hum))*0.005,
#                       inverted_quadratic(temp,-1.48e-01, 9.16, 37.73)+(1-(1.2248 + 0.2673*hum))*0.01))
# }

muv_th <- function(temp, c, T0, Tm, hum){
        return(ifelse(hum <= 1,inverted_quadratic(temp, c, T0, Tm)+(1-(0.01256 + 2.00893*hum))*0.005,
                      inverted_quadratic(temp, c, T0, Tm)+(1-(1.2248 + 0.2673*hum))*0.01))
}

# If humidity was not affecting
# muv_th <- function(temp, c, T0, Tm, hum){
#         
#                       inverted_quadratic(temp, c, T0, Tm)
# }


#' Larvae mortality rate
#' 
#' function of temperature and rain
#' @param t daily mean temperature
#' @param pr amount of rainfall summed over the prior two-week period
#' Note: anonymous function theta and temperature dependent mortality included
# mul_tr_t <- function(t,pr,pr_c=30,mu_al=0.001,mul_b0=2.32,mul_b1=-4.19e-1,mul_b2=2.73e-2,mul_b3=-7.53e-4,mul_b4=7.50e-6){
#         theta <- function(r,r_c=30){
#                 return(ifelse(r - r_c > 0 | r - r_c == 0, 1, 0))
#         }
#         mul_t <- function(x){
#                 return(poly(t, mul_b0, mul_b1, mul_b2, mul_b3, mul_b4))
#         }
#         # to stop rain impact remove theta(pr)
#         return(mul_t(t)*(1+mu_al*(pr-pr_c)*theta(pr)))
# }
#Now using quadratic fit
mul_tr_t <- function(t, c, T0, Tm, pr, pr_c=30,mu_al=0.001){
        theta <- function(r,r_c=30){
                return(ifelse(r - r_c > 0 | r - r_c == 0, 1, 0))
        }
        mul_t <- function(t, c, T0, Tm){
                return(inverted_quadratic(t, c, T0, Tm))
        }
        # to stop rain impact remove theta(pr)
        return(mul_t(t, c, T0, Tm)*(1+mu_al*(pr-pr_c)*theta(pr)))
}


#' Pupae mortality rate
#' 
#' Polynomial fit to empirical data
#' 
#' function of temperature
#' @param t daily mean temperature
#' @return return daily pupae mortality rate
#' Note: Liu Helmersson (2019) shows mup_b3 to be 4.39e-7
# mup_tr_t <- function(t,pr,pr_c=30,mu_ap=0.001,mup_b0=4.25e-1,mup_b1=-3.25e-2,mup_b2=7.06e-4,mup_b3=-4.39e-7){
#         theta <- function(r,r_c=30){
#                 return(ifelse(r - r_c > 0 | r - r_c == 0, 1, 0))
#         }
#         mup_t <- function(x) {
#                 return(poly(t, mup_b0, mup_b1, mup_b2, mup_b3))
#         }
#         # to stop rain impact remove theta(pr)
#         return(mup_t(t)*(1+mu_ap*(pr-pr_c)*theta(pr)))
# }
#Now using quadratic fit
mup_tr_t <- function(t, c, T0, Tm, pr, pr_c=30,mu_al=0.001){
        theta <- function(r,r_c=30){
                return(ifelse(r - r_c > 0 | r - r_c == 0, 1, 0))
        }
        mup_t <- function(t, c, T0, Tm){
                return(inverted_quadratic(t, c, T0, Tm))
        }
        # to stop rain impact remove theta(pr)
        return(mup_t(t, c, T0, Tm)*(1+mu_al*(pr-pr_c)*theta(pr)))
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
K_t <- function(pr,k0=5,k1=30,k2=0.1){
        return((k0*pr/(k1+pr))+k2)
}

#' The function fillGap is used to fill missing values in a climate variable dataset
#' 
#' @param clim_var a vector of climate variable values.
#' @param t_step a scalar representing the time step of the data.
#' @param dtr a vector of differences between consecutive values of the climate variable.
#' @param option a string that indicates the type of climate variable. It can be "temp", "hum", "rain", or "rain14".
# 
# fillGap <- function(clim_var,t_step,option){
#         if (option=="temp" | option=="hum" | option=="rain" | option=="rain14"){
#                 return(c(rep(clim_var[1:length(clim_var)-1],each=1/t_step),clim_var[length(clim_var)]))
#         }
# }

fillGap <- function(clim_var,t_step, dtr, option){
        if (option=="temp") {
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

        } else if (option=="hum" | option=="rain" | option=="rain14") {
                return(c(rep(clim_var[1:length(clim_var)-1],each=1/t_step),clim_var[length(clim_var)]))
        }
}

#' function to convert daily climate data to required timesteps and rearrange years for GCMs
#' 
#' @param climdata: a dataframe containing the climate data. It should contain the following 
#'                  columns: "gcm" (general circulation model), "rcp" (representative concentration pathway), 
#'                  "meanTemp" (mean temperature), "SVPD" (specific vapor pressure deficit), "prec" (precipitation), 
#'                  "prec_14day" (14 day accumulated precipitation), "date" (date of the observation).
#' @param source: a string indicating the source of the data, it can be "Observed" or a specific GCM (General Circulation Model).
#' @param scenario: a string indicating the scenario, it can be "Historical" or a specific RCP (Representative Concentration Pathway).
#' @param years: a vector containing the years of the data.
#' @param t_step: a scalar indicating the time step of the data.
#' @param pts: a scalar indicating the number of points of the data. 
#' 
#' 
climProcess <- function(climdata, source, scenario, years, t_step, pts){
        
        dat1 <- subset(climdata, climdata$gcm==source &
                               climdata$rcp=="Historical")
        dat2 <- subset(climdata, climdata$gcm==source &
                               climdata$rcp==scenario)
        # Return empty df when Observed with RCP or ISIMIP GCM with observed combination
        if (length(dat2$gcm)==0){
                return(NA)       
        }
        
        dat <- rbind(dat1, dat2)
        
        tempData <- list()
        tempClimData <- list()
        
        tempData$temp <- fillGap(dat$meanTemp, t_step, dat$dtr, option = "temp")
        # for no sinusoid
        # tempData$temp <- fillGap(dat$meanTemp, t_step, option = "temp")
        tempData$hum <- fillGap(dat$SVPD, t_step, option = "hum")
        tempData$rain <- fillGap(dat$prec, t_step, option = "rain")
        tempData$rain14 <- fillGap(dat$prec_14day, t_step, option = "rain14")
        # Observed data is not available after 2016
        if (source=="Observed"){
                opts <- seq(1, length(dat$date), by = t_step)
                tempData$tdate <- as_datetime(dat$date[1] + opts -1)
        } else {
                tempData$tdate <- as_datetime(dat$date[1] + pts -1)
        }
        tempData$year <- year(tempData$tdate)
        tempClimData <- bind_rows(tempClimData, tempData) %>%
                dplyr::select(tdate, everything())
        
        # replace missing values with previous observation
        tempClimData$temp <- zoo::na.locf(tempClimData$temp)
        tempClimData$hum <- zoo::na.locf(tempClimData$hum)
        tempClimData$rain <- zoo::na.locf(tempClimData$rain)
        tempClimData$rain14 <- zoo::na.locf(tempClimData$rain14)
        # no rearrangement for Observed data
        if (source == "Observed"){
                rearrangedData <- tempClimData
        } else {
                rearrangedData <- left_join(data.frame(year = years),    # Reorder data frame
                                            tempClimData,
                                            by = "year")
        }
        return(rearrangedData)
}
