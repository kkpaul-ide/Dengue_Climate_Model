## isimip_data_processing.R

# Script to read individual netCDF files containing climate data for 35 weather
# station locations and put all climate data into one dataframe

# Kishor Kumar Paul

# Set-up ------------------------------------------------------------------
# Restart R and source to main project directory

#check/install/load packages
mypackages<-c("tidyverse", "reshape2", "RcppRoll", "ncdf4", "zoo")

for (p in mypackages){
        if(!require(p, character.only = TRUE)){
                install.packages(p)
                library(p, character.only = TRUE)
        }
}

# Organize directories
basepath <- getwd()
isimipPath <- file.path(basepath, "data", "ISIMIP_data", "bangladesh")

# Create path for each model and scenario
gfdlRCP45_path <- file.path(isimipPath, "GFDL-ESM2M", "rcp45")
gfdlRCP85_path <- file.path(isimipPath, "GFDL-ESM2M", "rcp85")
hadgemRCP45_path <- file.path(isimipPath, "HadGEM2-ES", "rcp45")
hadgemRCP85_path <- file.path(isimipPath, "HadGEM2-ES", "rcp85")
ipslRCP45_path <- file.path(isimipPath, "IPSL-CM5A-LR", "rcp45")
ipslRCP85_path <- file.path(isimipPath, "IPSL-CM5A-LR", "rcp85")
mirocRCP45_path <- file.path(isimipPath, "MIROC5", "rcp45")
mirocRCP85_path <- file.path(isimipPath, "MIROC5", "rcp85")
gfdlHist_path <- file.path(isimipPath, "GFDL-ESM2M", "historical")
hadgemHist_path <- file.path(isimipPath, "HadGEM2-ES", "historical")
ipslHist_path <- file.path(isimipPath, "IPSL-CM5A-LR", "historical")
mirocHist_path <- file.path(isimipPath, "MIROC5", "historical")

# Path to observed data
obsPath <- "C:/Users/kisho/OneDrive - UNSW/PhD docs/Climate data/1975_2016"
obsFile <- "obs_climate_data.csv"

# Read weather station names 
df <- read.csv(file.path(basepath,"data","ListOf35BMDstations.csv"),
               header = TRUE, sep = ",")

# Make a list of divisions with weather stations
divs <- list(Mymensingh = "Mymensingh",
             Sylhet     = c("Sylhet","Srimangal"),
             Rajshahi   = c("Ishurdi", "Rajshahi", "Bogra"),
             Rangpur    = c("Syedpur", "Dinajpur", "Rangpur"),
             Dhaka      = c("Dhaka", "Faridpur", "Madaripur", "Tangail"),
             Barisal    = c("Barisal", "Bhola", "Khepupara", "Patuakhali"),
             Khulna     = c("Chuadanga", "Jessore", "Khulna", "Mongla", 
                            "Satkhira"),
             Chittagong = c("Chittagong (AP)", "Chandpur", 
                            "Chittagong (City)", "Comilla", "Cox's Bazar",
                            "Feni", "Hatiya", "Kutubdia", "Maijdee Court", 
                            "Rangamati", "Sandwip", "Sitakunda", "Teknaf"))


# Script parameters -------------------------------------------------------
# To be changed by user
reload <- FALSE
saveResults <- FALSE

# generic functions

#' kelvin to celsius
#' 
#' Convert temperature from kelvin to celsius
#' 
#'  @param temp temperature value in kelvin
#'  @return temperature value in celsius
#'   
KelvinCelsius <- function(temp){
        tempincelsius <- temp - 273.15
        return(tempincelsius)
}

#' listTodf
#' 
#' Convert a list to dataframe
#' 
#' @param lis the list to be converted
#' @return ...
#' 
listTodf <- function(lis){
        lis <- lis %>% # convert the list to a dataframe
                unlist(recursive = FALSE) %>% 
                enframe() 
}

#' isimipHist
#' 
#' function to extract temp, hum, prec data from three dimensional ISIMIP 
#' netCDF files and make two dimensional df for historical period
#' 
#' @param filepath File path to a ISIMIP model
#'  
#' @return returns df with mean temp, relative humidity and precepitation 
#' variables for a ISIMIP models

isimipHist <- function(filepath){
        # Read file names with daily mean temperature 
        files <- list.files(path = file.path(filepath, "tas"), 
                            pattern = "*_station.nc")
        station_meanTemp <- vector("list", length(files))
        for (i in 1:length(files)){
                # open nc files using ncdf4 package
                nc <- nc_open(file.path(filepath, "tas", files[i])) 
                # extract temperature variable data from a nc file
                station_meanTemp[[i]] <- ncvar_get(nc, "tasAdjust") 
        }
        # convert temperature unit using KelvinCelsius function
        station_meanTemp <- lapply(station_meanTemp, KelvinCelsius)
        
        # Read file names with daily mean temperature 
        files <- list.files(path = file.path(filepath, "tasmax"), 
                            pattern = "*_station.nc")
        station_maxTemp <- vector("list", length(files))
        for (i in 1:length(files)){
                # open nc files using ncdf4 package
                nc <- nc_open(file.path(filepath, "tasmax", files[i])) 
                # extract temperature variable data from a nc file
                station_maxTemp[[i]] <- ncvar_get(nc, "tasmaxAdjust") 
        }
        # convert temperature unit using KelvinCelsius function
        station_maxTemp <- lapply(station_maxTemp, KelvinCelsius)
        
        # Read file names with daily mean temperature 
        files <- list.files(path = file.path(filepath, "tasmin"), 
                            pattern = "*_station.nc")
        station_minTemp <- vector("list", length(files))
        for (i in 1:length(files)){
                # open nc files using ncdf4 package
                nc <- nc_open(file.path(filepath, "tasmin", files[i])) 
                # extract temperature variable data from a nc file
                station_minTemp[[i]] <- ncvar_get(nc, "tasminAdjust") 
        }
        # convert temperature unit using KelvinCelsius function
        station_minTemp <- lapply(station_minTemp, KelvinCelsius)
        
        # Read file names with daily relative humidity 
        files <- list.files(path = file.path(filepath, "hurs"), 
                            pattern = "*_station.nc")
        station_relHum <- vector("list", length(files))
        for (i in 1:length(files)){
                # open nc files using ncdf4 package
                nc <- nc_open(file.path(filepath, "hurs", files[i]))
                # extract humidity variable data from a nc file
                station_relHum[[i]] <- ncvar_get(nc, "hursAdjust")
        }
        
        # Read file names with daily precipitation 
        files <- list.files(path = file.path(filepath, "pr"), 
                            pattern = "*_station.nc")
        station_prec <- vector("list", length(files))
        for (i in 1:length(files)){
                # open nc files using ncdf4 package
                nc <- nc_open(file.path(filepath, "pr", files[i]))
                # extract precipitations variable data from a nc file
                station_prec[[i]] <- ncvar_get(nc, "prAdjust")
        }
        
        # calculate DTR from max and min temperature
        station_dtr <- mapply('-', station_maxTemp, station_minTemp, 
                              SIMPLIFY = FALSE)
        
        # use listTodf function convert list containing
        # temp data to a df
        station_meanTemp.df <- listTodf(station_meanTemp) %>% 
                rename(meanTemp = value) %>% 
                select(meanTemp)
        station_maxTemp.df <- listTodf(station_maxTemp) %>% 
                rename(maxTemp = value) %>% 
                select(maxTemp)
        station_minTemp.df <- listTodf(station_minTemp) %>% 
                rename(minTemp = value) %>% 
                select(minTemp)
        station_dtr.df <- listTodf(station_dtr) %>% 
                rename(dtr = value) %>% 
                select(dtr)
        station_relHum.df <- listTodf(station_relHum) %>% 
                rename(relHum = value) %>% 
                select(relHum)
        station_prec.df <- listTodf(station_prec) %>%
                rename(prec = value) %>%
                select(prec)
        
        
        # Column bind all temp and VC variables into one df
        df <- cbind(station_maxTemp.df, station_meanTemp.df, station_minTemp.df,
                    station_dtr.df, station_relHum.df,station_prec.df)
        return(df)
}

#' isimipRCP
#' 
#' function to extract temp, hum, prec data from three dimensional ISIMIP 
#' netCDF files and make two dimensional df for 2006-2099
#' 
#' @param filepath File path to a ISIMIP model
#'  
#' @return returns df with mean, max, min and DTR variables for a ISIMIP
#' model

isimipRCP <- function(filepath){
        # Read file names with daily mean temperature 
        files <- list.files(path = file.path(filepath, "tas"), 
                            pattern = "*_station.nc")
        station_meanTemp <- vector("list", length(files))
        for (i in 1:length(files)){
                # open nc files using ncdf4 package
                nc <- nc_open(file.path(filepath, "tas", files[i]))
                # extract temperature variable data from a nc file
                station_meanTemp[[i]] <- ncvar_get(nc, "tasAdjust") 
        }
        # Remove data beyond 31 Dec 2099 in some files
        station_meanTemp[[10]] <- station_meanTemp[[10]][,1:3287]
        # convert temperature unit using function from functions.R
        station_meanTemp <- lapply(station_meanTemp, KelvinCelsius) 
        
        # Read file names with daily max temperature 
        files <- list.files(path = file.path(filepath, "tasmax"), 
                            pattern = "*_station.nc")
        station_maxTemp <- vector("list", length(files))
        for (i in 1:length(files)){
                # open nc files using ncdf4 package
                nc <- nc_open(file.path(filepath, "tasmax", files[i]))
                # extract temperature variable data from a nc file
                station_maxTemp[[i]] <- ncvar_get(nc, "tasmaxAdjust")
        }
        # Remove data beyond 31 Dec 2099 in some files
        station_maxTemp[[10]] <- station_maxTemp[[10]][,1:3287]
        # convert temperature unit using function from functions.R
        station_maxTemp <- lapply(station_maxTemp, KelvinCelsius)
        
        # Read file names with daily min temperature 
        files <- list.files(path = file.path(filepath, "tasmin"), 
                            pattern = "*_station.nc")
        station_minTemp <- vector("list", length(files))
        for (i in 1:length(files)){
                # open nc files using ncdf4 package
                nc <- nc_open(file.path(filepath, "tasmin", files[i]))
                # extract temperature variable data from a nc file
                station_minTemp[[i]] <- ncvar_get(nc, "tasminAdjust")
        }
        # Remove data beyond 31 Dec 2099 in some files
        station_minTemp[[10]] <- station_minTemp[[10]][,1:3287]
        # convert temperature unit using function from functions.R
        station_minTemp <- lapply(station_minTemp, KelvinCelsius)
        
        # calculate DTR from max and min temperature
        station_dtr <- mapply('-', station_maxTemp, station_minTemp, 
                              SIMPLIFY = FALSE)
        
        # Read file names with daily relative humidity 
        files <- list.files(path = file.path(filepath, "hurs"), 
                            pattern = "*_station.nc")
        station_relHum <- vector("list", length(files))
        for (i in 1:length(files)){
                # open nc files using ncdf4 package
                nc <- nc_open(file.path(filepath, "hurs", files[i]))
                # extract temperature variable data from a nc file
                station_relHum[[i]] <- ncvar_get(nc, "hursAdjust")
        }
        # Remove data beyond 31 Dec 2099 in some files
        station_relHum[[10]] <- station_relHum[[10]][,1:3287]
        
        # Read file names with daily precipitation 
        files <- list.files(path = file.path(filepath, "pr"), 
                            pattern = "*_station.nc")
        station_prec <- vector("list", length(files))
        for (i in 1:length(files)){
                # open nc files using ncdf4 package
                nc <- nc_open(file.path(filepath, "pr", files[i]))
                # extract temperature variable data from a nc file
                station_prec[[i]] <- ncvar_get(nc, "prAdjust")
        }
        # Remove data beyond 31 Dec 2099 in some files
        station_prec[[10]] <- station_prec[[10]][,1:3287]
        
        # use listTodf function to convert list containing
        # temp data to a df
        station_meanTemp.df <- listTodf(station_meanTemp) %>% 
                rename(meanTemp = value) %>% 
                select(meanTemp)
        station_maxTemp.df <- listTodf(station_maxTemp) %>% 
                rename(maxTemp = value) %>% 
                select(maxTemp)
        station_minTemp.df <- listTodf(station_minTemp) %>% 
                rename(minTemp = value) %>% 
                select(minTemp)
        station_dtr.df <- listTodf(station_dtr) %>% 
                rename(dtr = value) %>% 
                select(dtr)
        station_relHum.df <- listTodf(station_relHum) %>% 
                rename(relHum = value) %>% 
                select(relHum)
        station_prec.df <- listTodf(station_prec) %>%
                rename(prec = value) %>%
                select(prec)
        
        # Column bind all temp and VC variables into one df
        df <- cbind(station_maxTemp.df, station_meanTemp.df, station_minTemp.df, 
                    station_dtr.df, station_relHum.df, station_prec.df)
        
        return(df)
}

# Extract data for ISIMIP models-------------------------------------------
if (reload) {
        load(file.path(basepath, "data", "climateData_station.RData"))
} else {
        # Execute function to extract climate data using 
        # function isimipRCP.R & isimipHist.R
        gfdlHist_station <- isimipHist(gfdlHist_path)
        hadgemHist_station <- isimipHist(hadgemHist_path)
        ipslHist_station <- isimipHist(ipslHist_path)
        mirocHist_station <- isimipHist(mirocHist_path)
        gfdlRCP45_station <- isimipRCP(gfdlRCP45_path)
        gfdlRCP85_station <- isimipRCP(gfdlRCP85_path)
        hadgemRCP45_station <- isimipRCP(hadgemRCP45_path)
        hadgemRCP85_station <- isimipRCP(hadgemRCP85_path)
        ipslRCP45_station <- isimipRCP(ipslRCP45_path)
        ipslRCP85_station <- isimipRCP(ipslRCP85_path)
        mirocRCP45_station <- isimipRCP(mirocRCP45_path)
        mirocRCP85_station <- isimipRCP(mirocRCP85_path)
        
        # Create df of stations corresponding to number of days 
        # available for each model and RCP scenario/historical period
        # RCP scenarios: Temp data for 34,333 days (Jan 2006 to Dec 2099) available  
        station.df <- data.frame(station = rep(df$Station, times = 34333))
                
        # Historical: Temp data for 20,089 days (Jan 1951 to Dec 2005) available   
        hist.station.df <- data.frame(station = rep(df$Station, times = 20089))
        
        # create df of dates corresponding to RCP scenarios 
        
        date.df <- data.frame(date = rep(seq(as.Date("2006/1/1"),
                                             as.Date("2099/12/31"), "days"), each=35))
        
        # create df of dates corresponding to historical period 
        
        hist.date.df <- data.frame(date = rep(seq(as.Date("1951/1/1"),
                                             as.Date("2005/12/31"), "days"), each=35))
        
        # column bind district, date and VC for each model
        
        gfdlRCP45_df <- cbind(station.df, date.df, gfdlRCP45_station) %>% 
                mutate(gcm = "GFDL-ESM2M", rcp = "RCP 4.5") 
        gfdlRCP85_df <- cbind(station.df, date.df, gfdlRCP85_station) %>% 
                mutate(gcm = "GFDL-ESM2M", rcp = "RCP 8.5")
        hadgemRCP45_df <- cbind(station.df, date.df, hadgemRCP45_station) %>% 
                mutate(gcm = "HadGEM2-ES", rcp = "RCP 4.5") 
        hadgemRCP85_df <- cbind(station.df, date.df, hadgemRCP85_station) %>% 
                mutate(gcm = "HadGEM2-ES", rcp = "RCP 8.5")
        ipslRCP45_df <- cbind(station.df, date.df, ipslRCP45_station) %>% 
                mutate(gcm = "IPSL-CM5A-LR", rcp = "RCP 4.5")
        ipslRCP85_df <- cbind(station.df, date.df, ipslRCP85_station) %>% 
                mutate(gcm = "IPSL-CM5A-LR", rcp = "RCP 8.5")
        mirocRCP45_df <- cbind(station.df, date.df, mirocRCP45_station) %>% 
                mutate(gcm = "MIROC5", rcp = "RCP 4.5")
        mirocRCP85_df <- cbind(station.df, date.df, mirocRCP85_station) %>% 
                mutate(gcm="MIROC5", rcp = "RCP 8.5")
        
        gfdlHist_df <- cbind(hist.station.df, hist.date.df, gfdlHist_station) %>% 
                mutate(gcm = "GFDL-ESM2M", rcp = "Historical") 
        hadgemHist_df <- cbind(hist.station.df, hist.date.df, hadgemHist_station) %>% 
                mutate(gcm="HadGEM2-ES", rcp = "Historical")
        ipslHist_df <- cbind(hist.station.df, hist.date.df, ipslHist_station) %>% 
                mutate(gcm="IPSL-CM5A-LR", rcp = "Historical")
        mirocHist_df <- cbind(hist.station.df, hist.date.df, mirocHist_station) %>% 
                mutate(gcm="MIROC5", rcp = "Historical")
        
        # read and process observed data
        obs_clim <- read.csv(file.path(obsPath, obsFile))
        obs <- obs_clim %>%
                mutate(date = as.Date(with(obs_clim,paste(year,month,day,sep="-")),"%Y-%m-%d")) %>%
                mutate(dtr = max_temp - min_temp) %>%
                mutate(avg_temp = na.locf(avg_temp)) %>% # replace NAs with previous temp
                select(-c(year, month, day)) %>%
                rename(meanTemp = avg_temp,
                       minTemp = min_temp,
                       maxTemp = max_temp,
                       relHum = avg_hum,
                       prec = tot_rain) %>%
                mutate(gcm = "Observed") %>%
                mutate(rcp = "Observed") %>%
                mutate(prec = prec/86400)
        
        # Remove intermediate dfs to free up space
        rm(date.df,hist.date.df,station.df,hist.station.df,mirocHist_station,ipslHist_station,hadgemHist_station,
           gfdlHist_station,mirocRCP85_station,mirocRCP45_station,ipslRCP85_station,ipslRCP45_station,
           gfdlRCP45_station,gfdlRCP85_station,hadgemRCP45_station,hadgemRCP85_station)
        
        # Bind data for each model into one df

        hist.station.df <- reduce(list(gfdlRCP45_df, gfdlRCP85_df, hadgemRCP45_df, hadgemRCP85_df, ipslRCP45_df, 
                                       ipslRCP85_df, mirocRCP45_df, mirocRCP85_df, gfdlHist_df, hadgemHist_df, 
                                       ipslHist_df, mirocHist_df, obs), full_join)
        
        # Remove intermediate dfs to free up space
        rm (ipslRCP45_df, ipslRCP85_df, gfdlRCP45_df, gfdlRCP85_df, hadgemRCP45_df, hadgemRCP85_df,
            mirocRCP45_df, mirocRCP85_df,gfdlHist_df,hadgemHist_df,ipslHist_df,mirocHist_df)
        
        # Take rolling sum of precipitation in previous 14 days and convert unit of precipitation from 
        # kg m-2 s-1 to mm/day (1 kg m-2 s-1 = 86400 mm/day)
        clim.station.df <- hist.station.df %>%
                mutate(prec = prec*86400) %>%
                arrange(gcm, rcp, station, date) %>%
                group_by(gcm, rcp, station) %>%
                mutate(prec_14day = roll_sum(prec, 14, align = "right", fill = NA)) 
                
        
        # calculate SVPD as in Caldwell et. al. 2021 & http://cronklab.wikidot.com/calculation-of-vapour-pressure-deficit
        clim.station.df$SVPD <- (1-(clim.station.df$relHum/100))*610.7*10^(7.5*clim.station.df$meanTemp/(237.3+clim.station.df$meanTemp))/1000
                
        
        # Change variable positions for easy readability
        clim.station.df <- clim.station.df[ , c("gcm", "rcp", "station", "date", "meanTemp", "maxTemp",
                                                "minTemp", "dtr", "relHum", "prec", "prec_14day", "SVPD")]


        if (saveResults) {
                save(clim.station.df, file=file.path(basepath, "data",  
                                                     "climateData_station.RData"))
        }
}





