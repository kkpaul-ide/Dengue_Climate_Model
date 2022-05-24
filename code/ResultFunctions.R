#' function to aggregate model simulation output and determine seroprev on a specific date  


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
                summarise(infections = sum(newInfections))

        monthInfections <- dat %>%
                group_by(year, month) %>%
                summarise(infections = sum(newInfections))

        # Sero prevalence end of 2013
        seroprev <- dat$seroPrev[pg$pts[pg$pts_times == as_datetime(date)]]
        
        return(list(dat = dat, annInfections = annInfections,
                    monthInfections = monthInfections, seroprev = seroprev))
        
}


gcmPlotMonthly <- function (dat, gcm) {
        ggplot(dat %>% 
                       filter(dat$Model==gcm) %>% 
                       group_by(year, month) %>% 
                       summarise(Infections = mean(infections), 
                                 minInfections = min(infections), 
                                 maxInfections = max(infections)),
               aes(x = year + month/12)) + 
                geom_ribbon(aes(ymin = minInfections, ymax = maxInfections), fill = "grey") +
                geom_line(aes(y = Infections),color = "red") +
                scale_x_continuous(breaks = seq(pg$startYear, pg$endYear),
                                   labels = EveryNth(seq(pg$startYear, pg$endYear), 5, inverse = TRUE),
                                   limits = c(pg$startYear, pg$endYear)) +
                xlab(paste(gcm)) + ylab("Dengue infections per month") + 
                coord_cartesian(xlim = c(pg$startYear, pg$endYear), ylim = c(0, 6e5)) +
                PlotOptions()
}