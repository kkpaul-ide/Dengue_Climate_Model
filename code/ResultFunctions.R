#' Function to aggregate model simulation output and determine seroprev on a specific date.
#' The AdjustResults() function in R takes three inputs:
#' 
#' @param dat a data frame containing the simulation results
#' @param date a date string in the format of yyyy-mm-dd
#' @param pg a list containing the parameters for the simulation

AdjustResults <- function(dat, date, pg) {
        
        if (pg$model[1] == "base"){
                # dat$seroPrev <- (dat$Nh-dat$s_h)*100/dat$Nh
                dat$seroPrev <- (dat$r_h1+dat$r_h2+dat$r)*100/dat$Nh
        } else if (pg$model[1] == "caldwell") {
                dat$seroPrev <- (dat$e_h + dat$i_h + dat$r)*100/dat$Nh
        }
        
        dat <- aggregate(dat, list(rep(1:(nrow(dat) %/% (1/pg$dt) + 1), each = 1/pg$dt, 
                                       len = nrow(dat))), mean)[-1]
        dat$newInfections <- (1/pg$dt)*dat$newInfections
        dat$newInfections1 <- (1/pg$dt)*dat$newInfections1
        dat$newInfections2 <- (1/pg$dt)*dat$newInfections2
        dat$days <- 1:nrow(dat)
        dat$date <- pg$start_date + dat$days-1
        dat$month <- month(dat$date)
        dat$year <- year(dat$date)

        annInfections <- dat %>%
                group_by(year) %>%
                summarise(infections = sum(newInfections), na.rm = TRUE)

        monthInfections <- dat %>%
                group_by(year, month) %>%
                summarise(infections = sum(newInfections))

        # Sero prevalence end of 2013
        seroprev <- dat$seroPrev[pg$pts[pg$pts_times == as_datetime(date)]]
        
        return(list(dat = dat, annInfections = annInfections,
                    monthInfections = monthInfections, seroprev = seroprev))
        
}

