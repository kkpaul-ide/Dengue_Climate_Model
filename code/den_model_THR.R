#' Two strain differential dengue transmission model
#' 
#' Model function
#' 
#' @param t Integer value specifying the number of days the model is 
#' run for.
#' @param state Named vector which contains the initial values
#' of each of the mosquito and human compartment. 
#' @param params Named vector/dataframe of numeric values which contains 
#' the time independent input parameters for the model. 
#' @return Return a dataframe containing values of each compartment
#' at specified time points

den_model_THR <- function(t, state, params) {
        with(as.list(c(state,params)), {
                Nv <- Sv + Ev1 + Ev2 + Iv1 + Iv2
                Nl <- Sl + Il1 + Il2
                Nh <- Sh + Eh1 + Eh2 + Ih1 + Ih2 + CIh1 + CIh2 + Rh1 + Rh2 + 
                        Eh12 + Eh21 + Ih12 + Ih21 + Rh
               
                 # ODEs for  mosquito populations 
                dSl <- (muv_th(temp[t], hum[t], timestep)*Nv + mul_r(rain[t])*Nl)*(1 - VTR*(Iv1 + Iv2)/Nv) -
                        (PMR(temp[t]) + PMR(temp[t]))*Sl -  mul_r(rain[t])*Sl 
                dIl1 <- (muv_th(temp[t], hum[t], timestep)*Nv + mul_r(rain[t])*Nl) * VTR*Iv1/Nv - PMR(temp[t])*Il1 - 
                        mul_r(rain[t])*Il1
                dIl2 <- (muv_th(temp[t], hum[t], timestep)*Nv + mul_r(rain[t])*Nl) * VTR*Iv2/Nv - PMR(temp[t])*Il2 - 
                        mul_r(rain[t])*Il2
                dSv <- (PMR(temp[t]) + PMR(temp[t]))*Sl - a(temp[t])*pMI(temp[t])*(Ih1 + Ih2 + ADE*Ih12 + ADE*Ih21)/Nh * Sv - 
                        muv_th(temp[t], hum[t], timestep)*Sv
                dEv1 <- a(temp[t]) * pMI(temp[t]) * (Ih1 + ADE*Ih21)/Nh*Sv - VIR(temp[t])*Ev1 - muv_th(temp[t], hum[t], timestep)*Ev1
                dEv2 <- a(temp[t]) * pMI(temp[t]) * (Ih2 + ADE*Ih12)/Nh*Sv - VIR(temp[t])*Ev2 - muv_th(temp[t], hum[t], timestep)*Ev2
                dIv1 <- PMR(temp[t])*Il1 + VIR(temp[t])*Ev1 - muv_th(temp[t], hum[t], timestep)*Iv1
                dIv2 <- PMR(temp[t])*Il2 + VIR(temp[t])*Ev2 - muv_th(temp[t], hum[t], timestep)*Iv2
                
                # ODEs for human populations        
                dSh <- BR*(Nh/1000)/360 - a(temp[t])*pIM(temp[t])*(Iv1 + Iv2)/Nh * Sh - DR*(Sh/1000)/360
                dEh1 <- a(temp[t])*pIM(temp[t])*Iv1/Nh * Sh - (1.0/5.9)*Eh1 - DR*(Eh1/1000)/360
                dEh2 <- a(temp[t])*pIM(temp[t])*Iv2/Nh * Sh - (1.0/5.9)*Eh2 - DR*(Eh2/1000)/360
                dIh1 <- (1.0/5.9)*Eh1 - (1.0/5.0)*Ih1 - DR*(Ih1/1000)/360
                dIh2 <- (1.0/5.9)*Eh2 - (1.0/5.0)*Ih2 - DR*(Ih2/1000)/360
                dCIh1 <- (1.0/5.0)*Ih1 - (1.0/dCI)*CIh1 - DR*(CIh1/1000)/360
                dCIh2 <- (1.0/5.0)*Ih2 - (1.0/dCI)*CIh2 - DR*(CIh2/1000)/360
                dRh1 <- (1.0/dCI)*CIh1 - a(temp[t])*pIM(temp[t])*Iv2/Nh*Rh1 - DR*(Rh1/1000)/360
                dRh2 <- (1.0/dCI)*CIh2 - a(temp[t])*pIM(temp[t])*Iv1/Nh*Rh2 - DR*(Rh2/1000)/360
                dEh12 <- a(temp[t])*pIM(temp[t])*Iv2/Nh*Rh1 - (1.0/5.9)*Eh12 - DR*(Eh12/1000)/360
                dEh21 <- a(temp[t])*pIM(temp[t])*Iv1/Nh*Rh2 - (1.0/5.9)*Eh21 - DR*(Eh21/1000)/360
                dIh12 <- (1.0/5.9)*Eh12 - (1.0/5.0 + DIM)*Ih12 - DR*(Ih12/1000)/360
                dIh21 <- (1.0/5.9)*Eh21 - (1.0/5.0 + DIM)*Ih21 - DR*(Ih21/1000)/360      
                dRh <-    (1.0/5.0)*Ih12 + (1.0/5.0)*Ih21 - DR*(Rh/1000)/360 
                
                
                list(c(dSl, dIl1, dIl2, dSv, dEv1, dEv2, dIv1, dIv2, dSh, dEh1, 
                       dEh2, dIh1, dIh2, dCIh1, dCIh2, dRh1, dRh2, dEh12, dEh21,
                       dIh12, dIh21, dRh))
        })
}



# This is the general function for the Briere fit.
briere <- function(x, c, T0, Tm){
        if((x < T0) | (x > Tm))
                0.0
        else
                c*x*(x-T0)*sqrt(Tm-x)
}


# This is the general function for the inverted quadratic fit.
inverted_quadratic <- function(x, c, T0, Tm, timestep){
        if((x < T0) | (x > Tm))
                1.0/timestep
        else
                1.0/(c*(x-T0)*(x-Tm))
}


# adult mosquito mortality rate (1/adult lifespan)
# depend on temperature and humidity
muv_th <- function(temp, hum, timestep){
        if (hum <= 1){
                inverted_quadratic(temp,-1.48e-01, 9.16, 37.73, timestep)+(1-(0.01256 + 2.00893*hum))*0.005
        } else {
                inverted_quadratic(temp,-1.48e-01, 9.16, 37.73, timestep)+(1-(1.2248 + 0.2673*hum))*0.01
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
