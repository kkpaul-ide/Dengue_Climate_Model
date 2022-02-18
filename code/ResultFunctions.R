#' function to aggregate model simulation output and determine seroprev on a specific date  


AdjustResults <- function(dat, source, scenario, pg) {
        
        if (pg$model[1] == "base"){
                dat$seroPrev <- (dat$Nh-dat$s_h)*100/dat$Nh
        } else if (pg$model[1] == "caldwell") {
                dat$seroPrev <- (dat$e_h + dat$i_h + dat$r)*100/dat$Nh
        }
        
        dat <- aggregate(dat, list(rep(1:(nrow(dat) %/% (1/pg$dt) + 1), each = 1/pg$dt, 
                                       len = nrow(dat))), mean)[-1]
        dat$newInfections <- (1/pg$dt)*dat$newInfections
        dat$days <- 1:nrow(dat)
        dat$date <- pg$start_date + dat$days-1
        dat$month <- month(dat$date)
        dat$year <- year(dat$date)
        
        
        annInfections <- aggregate(dat$newInfections, 
                                   list(rep(1:(nrow(dat) %/% 365 + 1), each = 365, 
                                            len = nrow(dat))), sum)[-1] %>%
                mutate(year = pg$startYear+1:nrow(.)) %>%
                rename(infections = x)
        
        annInfections <- dat %>%
                group_by(year) %>%
                summarise(infections = sum(newInfections))
        
        monthInfections <- dat %>% 
                group_by(year, month) %>%
                summarise(infections = sum(newInfections))
        
        # Sero prevalence end of 2103
        seroprev <- dat$seroPrev[pg$pts[pg$pts_times == as_datetime("2013-12-31")]]
        
        return(list(dat = dat, annInfections = annInfections, 
                    monthInfections = monthInfections, seroprev = seroprev))
        
}
