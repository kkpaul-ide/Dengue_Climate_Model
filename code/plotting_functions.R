
#' The extractorRData function is used to extract an object from an R data file
#' @param file: a string specifying the file path of the R data file
#' @param object: a string specifying the name of the object to be extracted from the file

extractorRData <- function(file, object){
        E <- new.env()
        load(file = file, envir = E)
        return(get(object, envir = E, inherits = F))
}

#' plot_MonInf is an R function that creates a plot of monthly infections for a given 
#' baseline seroepidemiological scenario, and for two specified time periods and two 
#' representative concentration pathways (RCPs). The function takes four input arguments:
#' 
#' 
#' @param baseline_sero a character string that specifies the baseline seroepidemiological scenario to be plotted
#' @param period1 a character string that specifies the first time period to be plotted (default is "2030_2049")
#' @param period2 a character string that specifies the second time period to be plotted (default is "2080_2099")
#' @param RCP1 a character string that specifies the first RCP to be plotted (default is "RCP 4.5")
#' @param RCP2 a character string that specifies the second RCP to be plotted (default is "RCP 8.5")
#' 
#' The function creates a list of plots, where each plot shows the median monthly infections for a given GCM, RCP 
#' and time period. These plots are then combined into a single plot using the plot_grid function. The final plot 
#' is a grid of four subplots arranged vertically, where each subplot corresponds to one GCM and RCP combination. 
#' The final plot also includes a legend indicating the GCM used for each plot. It is important to note that the 
#' function uses the extractorRData function to load the data from the Rda file sim_summary.Rda in the folder 
#' specified by the path simFolder, which is constructed by the function using the input arguments. 
#' The function also requires the ggplot2 library, the scales library and the stringr library.

plot_MonInf <- function(baseline_sero, period1 = "2030_2049", period2 = "2080_2099", period3 = "1995_2014",
                        RCP1 = "RCP 4.5", RCP2 = "RCP 8.5", RCP3 = "Observed"){
        
        periods <- c(period1, period2)
        RCPs <- c(RCP1, RCP2)
        tempdat <- list()
        
        for (i in 1:length(periods)){
                for (j in 1:length(RCPs)){
                        
                        simFolder <- file.path(projectFolder, "savedsims", paste0(periods[i]), paste0(RCPs[j]), paste0(baseline_sero))
                        monthInfectionsDf <- extractorRData(file.path(simFolder,"sim_summary.Rda"), 'monthInfectionsDf')
                        monthInfectionsDf[c('Sim', 'GCM', 'RCP')] <- str_split_fixed(monthInfectionsDf$model, ' ', 3)
                        name <- paste(periods[i],RCPs[j])
                        tempdat[[name]] <- monthInfectionsDf %>% group_by(year, month, GCM) %>%
                                summarise(medNewInf = median(infections)) %>%
                                # replace really small values with 1
                                mutate(medNewInf = ifelse(medNewInf<1, 1, medNewInf))
                }
        }
        dat <- unnest(enframe(tempdat))
        ## replacing the year values only for plotting purpose
        dat$year <- rep(1995:2014, times=4, each=48)
        
        simFolder <- file.path(projectFolder, "savedsims", paste0(period3), paste0(RCP3))
        monthInfectionsDf <- extractorRData(file.path(simFolder,"sim_summary.Rda"), 'monthInfectionsDf')
        monthInfectionsDf[c('Sim', 'GCM', 'RCP')] <- str_split_fixed(monthInfectionsDf$model, ' ', 3)
        dat1 <- monthInfectionsDf %>% group_by(year, month, GCM) %>%
                summarise(medNewInf = median(infections)) %>%
                # replace really small values with 1
                mutate(medNewInf = ifelse(medNewInf<1, 1, medNewInf))

        yearPoints <- dat1$year+dat1$month/12
        ptsPlots <- yearPoints[seq(6, length(yearPoints), 12)]
        yearLabels <- paste0(1:20)
        
        plot1 <- ggplot(dat) +
                geom_line(aes(x = year+month/12, y=medNewInf, color = GCM))+
                xlab("Year") + ylab("Monthly infections") +
                facet_wrap(~name, ncol = 1, scales = "free_x")+
                scale_x_continuous(breaks = ptsPlots,
                                   labels = yearLabels) +
                geom_line(dat1, mapping = aes(x = year+month/12, y=medNewInf))+
                theme(legend.title = element_blank(), legend.position="none",
                      text = element_text(face = "bold",size=18,colour="black"),
                      panel.background = element_blank(),
                      axis.line = element_line(colour = "black"))
        
        plot2 <- ggplot(dat) +
                geom_line(aes(x = year+month/12, y=medNewInf, color = GCM))+
                xlab("Year") + ylab("Monthly infections (log)") +
                scale_y_log10()+
                facet_wrap(~name, ncol = 1, scales = "free_x")+
                scale_x_continuous(breaks = ptsPlots,
                                   labels = yearLabels) +
                geom_line(dat1, mapping = aes(x = year+month/12, y=medNewInf))+
                theme(legend.title = element_blank(), legend.position="none",
                        text = element_text(face = "bold",size=18,colour="black"),
                      panel.background = element_blank(),
                      axis.line = element_line(colour = "black"))
        
        legend <- get_legend(plot1 + guides(color = guide_legend(nrow = 1)) + theme(legend.position = "top"))
        
        plot3 <- plot_grid(plot1, plot2, labels = c("A", "B"),
                          hjust = -4, label_size=18)
        
        plot <- plot_grid(legend, plot3, ncol = 1, rel_heights = c(.1, 1))
        
        return(plot)
}


#' The plot_SeroPrev function takes in several inputs and generates a plot that compares the seroprevalence 
#' (the proportion of a population that has been infected with a disease and has developed antibodies against it) 
#' over time in two different periods and under two different RCP scenarios. 
#' 
#' @param baseline_sero a string that specifies the baseline seroprevalence used in the simulation
#' @param period1 a string that specifies the first time period of the simulation
#' @param period2 a string that specifies the second time period of the simulation
#' @param RCP1 a string that specifies the first RCP scenario used in the simulation
#' @param RCP2 a string that specifies the second RCP scenario used in the simulation
#' @param year_range a string that specifies the year range of the observations.#' 
#' 

plot_SeroPrev <- function(baseline_sero, period1 = "2030_2049", period2 = "2080_2099", period3 = "1995_2014",
                          RCP1 = "RCP 4.5", RCP2 = "RCP 8.5", RCP3 = "Observed"){
        
        periods <- c(period1, period2)
        RCPs <- c(RCP1, RCP2)
        tempdat <- list()
        
        for (i in 1:length(periods)){
                for (j in 1:length(RCPs)){
                        
                        simFolder <- file.path(projectFolder, "savedsims", paste0(periods[i]), paste0(RCPs[j]), paste0(baseline_sero))
                        datDf <- extractorRData(file.path(simFolder,"sim_summary.Rda"), 'datDf')
                        # Create Model and Set variable from combined model variable
                        datDf[c('Sim', 'GCM', 'RCP')] <- str_split_fixed(datDf$model, ' ', 3)
                        name <- paste(periods[i],RCPs[j])
                        tempdat[[name]] <- datDf %>% group_by(days, GCM) %>%
                                summarise(medPrev = median(seroPrev))
                        rm(datDf)
                }
        }
        dat <- unnest(enframe(tempdat))
        
        simFolder <- file.path(projectFolder, "savedsims", paste0(period3), paste0(RCP3))
        datDf <- extractorRData(file.path(simFolder,"sim_summary.Rda"), 'datDf')
        # Create Model and Set variable from combined model variable
        datDf[c('Sim', 'GCM', 'RCP')] <- str_split_fixed(datDf$model, ' ', 3)
        
        dat1 <- datDf %>% group_by(days) %>%
                summarise(medPrev = median(seroPrev))
        
        
        # Plot specifications
        yearPoints <- seq(from = 1, to = 40, by = 2)
        ptsPlots <- yearPoints*182
        yearLabels <- paste(c(1:20))
  

        
        plot <- ggplot(dat) +
                geom_line(aes(x = days, y=medPrev, color = GCM))+
                xlab("") + ylab("Seroprevalence") +
                scale_y_continuous(limits = c(0,85))+
                scale_x_continuous(breaks = ptsPlots,
                                   labels = yearLabels) +
                facet_wrap(~name, ncol = 2, scales = "free_x")+
                geom_line(dat1, mapping = aes(days, medPrev))+
                theme(legend.title = element_blank(), legend.position="top",
                      text = element_text(face = "bold",size=18,colour="black"),
                      panel.background = element_blank(),
                      axis.line = element_line(colour = "black"))
        
        
}

#'The plot_annPeak function is used to plot boxplots of the annual peak number of infections 
#'in a simulation. It takes in the following parameters:
#'
#' @param baseline_sero a string that specifies the baseline seroprevalence used in the simulation
#' @param period1 a string that specifies the first time period of the simulation
#' @param period2 a string that specifies the second time period of the simulation
#' @param RCP1 a string that specifies the first RCP scenario used in the simulation
#' @param RCP2 a string that specifies the second RCP scenario used in the simulation
#' @param year_range a string that specifies the year range of the observations.
#'

plot_annPeak <- function(baseline_sero, period1 = "2030_2049", period2 = "2080_2099", 
                         RCP1 = "RCP 4.5", RCP2 = "RCP 8.5", year_range="1995_2014"){
        
        simFolder <- file.path(projectFolder, "savedsims", paste0(year_range), "Observed")
        # load summary file with output from 1000 sims
        datDf <- extractorRData(file.path(simFolder,"sim_summary.Rda"), 'datDf')
        # Create Model and Set variable from combined model variable
        datDf[c('Sim', 'GCM', 'RCP')] <- str_split_fixed(datDf$model, ' ', 3)
        df1 <-  aggregate(datDf$newInfections, list(datDf$year,datDf$Sim, datDf$GCM, datDf$RCP), FUN=max)
        rm(datDf)
        
        simFolder <- file.path(projectFolder, "savedsims", paste0(year_range), RCP1)
        # load summary file with output from 1000 sims
        datDf <- extractorRData(file.path(simFolder,"sim_summary.Rda"), 'datDf')
        # Create Model and Set variable from combined model variable
        datDf[c('Sim', 'GCM', 'RCP')] <- str_split_fixed(datDf$model, ' ', 3)
        df2 <-  aggregate(datDf$newInfections, list(datDf$year,datDf$Sim, datDf$GCM, datDf$RCP), FUN=max)  %>%
                mutate(Group.3 = paste0(year_range,"_",Group.3))
        rm(datDf)
        
        simFolder <- file.path(projectFolder, "savedsims", paste0(period1), paste0(RCP1), paste0(baseline_sero))
        datDf <- extractorRData(file.path(simFolder,"sim_summary.Rda"), 'datDf')
        # Create Model and Set variable from combined model variable
        datDf[c('Sim', 'GCM', 'RCP')] <- str_split_fixed(datDf$model, ' ', 3)
        df3 <- aggregate(datDf$newInfections, list(datDf$year,datDf$Sim, datDf$GCM, datDf$RCP), FUN=max) %>%
                mutate(Group.3 = paste0(period1,"_",Group.3))
        rm(datDf)
        
        simFolder <- file.path(projectFolder, "savedsims", paste0(period2), paste0(RCP1), paste0(baseline_sero))
        datDf <- extractorRData(file.path(simFolder,"sim_summary.Rda"), 'datDf')
        # Create Model and Set variable from combined model variable
        datDf[c('Sim', 'GCM', 'RCP')] <- str_split_fixed(datDf$model, ' ', 3)
        df4 <- aggregate(datDf$newInfections, list(datDf$year,datDf$Sim, datDf$GCM, datDf$RCP), FUN=max) %>%
                mutate(Group.3 = paste0(period2,"_",Group.3))
        rm(datDf)
        
        df <- rbind(df1, df2, df3, df4)
        tick_label <- c('GFDL-ESM2M', 'HadGEM2-ES', 'IPSL-CM5A-LR', 'MIROC5')
        # reverse the elements in levels to have the expected order of the bars in final plot
        df$Group.3 <- factor(df$Group.3, levels = c(as.character(unique(df$Group.3))))
        cols <- c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF")
        box_cols <- c("gray", cols, cols, cols)
        plot1 <- ggplot(df, aes(x = Group.3, y = x, fill = Group.3)) +
                geom_boxplot(alpha=0.3, outlier.shape = NA)+
                xlab("") + ylab("Annual peak infections") +
                scale_y_continuous(limits = c(0,2.5e5), labels = scales::scientific)+
                scale_x_discrete(breaks=c(as.character(unique(df$Group.3))),
                                 labels=c("Observed", tick_label, tick_label, tick_label))+
                geom_vline(xintercept = c(1.5,5.5,9.5), linetype = "longdash")+
                annotate("text", x = c(3.5,7.5,11.5), y = 2.5e5, 
                         label = c("1995_2014","2030-2049", "2080-2099"), size = 8)+
                scale_fill_manual(values=box_cols)+
                theme(text = element_text(face = "bold",size=18,colour="black"),
                      axis.text.x = element_text(angle = 45, hjust = 1),
                      panel.background = element_blank(),
                      axis.line = element_line(colour = "black"),
                      legend.position="none")
        
        rm(df, df2, df3, df4)
        
        simFolder <- file.path(projectFolder, "savedsims", paste0(year_range), RCP2)
        # load summary file with output from 1000 sims
        datDf <- extractorRData(file.path(simFolder,"sim_summary.Rda"), 'datDf')
        # Create Model and Set variable from combined model variable
        datDf[c('Sim', 'GCM', 'RCP')] <- str_split_fixed(datDf$model, ' ', 3)
        df2 <-  aggregate(datDf$newInfections, list(datDf$year,datDf$Sim, datDf$GCM, datDf$RCP), FUN=max)  %>%
                mutate(Group.3 = paste0(year_range,"_",Group.3))
        rm(datDf)
        
        simFolder <- file.path(projectFolder, "savedsims", paste0(period1), paste0(RCP2), paste0(baseline_sero))
        datDf <- extractorRData(file.path(simFolder,"sim_summary.Rda"), 'datDf')
        # Create Model and Set variable from combined model variable
        datDf[c('Sim', 'GCM', 'RCP')] <- str_split_fixed(datDf$model, ' ', 3)
        df3 <- aggregate(datDf$newInfections, list(datDf$year,datDf$Sim, datDf$GCM, datDf$RCP), FUN=max) %>%
                mutate(Group.3 = paste0(period1,"_",Group.3))
        rm(datDf)
        
        simFolder <- file.path(projectFolder, "savedsims", paste0(period2), paste0(RCP2), paste0(baseline_sero))
        datDf <- extractorRData(file.path(simFolder,"sim_summary.Rda"), 'datDf')
        # Create Model and Set variable from combined model variable
        datDf[c('Sim', 'GCM', 'RCP')] <- str_split_fixed(datDf$model, ' ', 3)
        df4 <- aggregate(datDf$newInfections, list(datDf$year,datDf$Sim, datDf$GCM, datDf$RCP), FUN=max) %>%
                mutate(Group.3 = paste0(period2,"_",Group.3))
        rm(datDf)
        
        df <- rbind(df1, df2, df3, df4)
        # reverse the elements in levels to have the expected order of the bars in final plot
        df$Group.3 <- factor(df$Group.3, levels = c(as.character(unique(df$Group.3))))
        plot2 <- ggplot(df %>% filter(Group.4 != RCP1), aes(x = Group.3, y = x, fill = Group.3)) +
                geom_boxplot(alpha=0.3, outlier.shape = NA)+
                xlab("") +
                scale_y_continuous(limits = c(0,2.5e5), labels = scales::scientific)+
                scale_x_discrete(breaks=c(as.character(unique(df$Group.3))),
                                 labels=c("Observed", tick_label, tick_label, tick_label))+
                geom_vline(xintercept = c(1.5,5.5, 9.5), linetype = "longdash")+
                annotate("text", x = c(3.5,7.5,11.5), y = 2.5e5, 
                         label = c("1995_2014","2030-2049", "2080-2099"), size = 8)+
                scale_fill_manual(values=box_cols)+
                theme(text = element_text(face = "bold",size=18,colour="black"),
                      axis.text.x = element_text(angle = 45, hjust = 1),
                      panel.background = element_blank(),
                      axis.line = element_line(colour = "black"),
                      axis.title.y =element_blank(),
                      legend.position="none")
        plot <- plot_grid(plot1, plot2, labels = c("A", "B"),
                           label_size=14)
        return(plot)
}

#'The plot_annInf function is used to generate box plots of the annual total infections under different climate scenarios.
#'
#' @param baseline_sero a string that specifies the baseline seroprevalence used in the simulation
#' @param period1 a string that specifies the first time period of the simulation
#' @param period2 a string that specifies the second time period of the simulation
#' @param RCP1 a string that specifies the first RCP scenario used in the simulation
#' @param RCP2 a string that specifies the second RCP scenario used in the simulation
#' @param year_range a string that specifies the year range of the observations.
#' 
plot_annInf <- function(baseline_sero, period1 = "2030_2049", period2 = "2080_2099", 
                         RCP1 = "RCP 4.5", RCP2 = "RCP 8.5", year_range="1995_2014"){
        
        simFolder <- file.path(projectFolder, "savedsims", paste0(year_range), "Observed")
        # load summary file with output from 1000 sims
        annInfectionsDf <- extractorRData(file.path(simFolder,"sim_summary.Rda"), 'annInfectionsDf')
        # Create Model and Set variable from combined model variable
        annInfectionsDf[c('Sim', 'GCM', 'RCP')] <- str_split_fixed(annInfectionsDf$model, ' ', 3)
        df1 <- annInfectionsDf
       
        simFolder <- file.path(projectFolder, "savedsims", paste0(year_range), RCP1)
        # load summary file with output from 1000 sims
        annInfectionsDf <- extractorRData(file.path(simFolder,"sim_summary.Rda"), 'annInfectionsDf')
        # Create Model and Set variable from combined model variable
        annInfectionsDf[c('Sim', 'GCM', 'RCP')] <- str_split_fixed(annInfectionsDf$model, ' ', 3)
        df2 <-  annInfectionsDf %>% mutate(GCM = paste0(year_range,"_",GCM))
        
        simFolder <- file.path(projectFolder, "savedsims", paste0(period1), paste0(RCP1), paste0(baseline_sero))
        # load summary file with output from 1000 sims
        annInfectionsDf <- extractorRData(file.path(simFolder,"sim_summary.Rda"), 'annInfectionsDf')
        # Create Model and Set variable from combined model variable
        annInfectionsDf[c('Sim', 'GCM', 'RCP')] <- str_split_fixed(annInfectionsDf$model, ' ', 3)
        df3 <- annInfectionsDf %>% mutate(GCM = paste0(period1,"_",GCM))
        
        simFolder <- file.path(projectFolder, "savedsims", paste0(period2), paste0(RCP1), paste0(baseline_sero))
        # load summary file with output from 1000 sims
        annInfectionsDf <- extractorRData(file.path(simFolder,"sim_summary.Rda"), 'annInfectionsDf')
        # Create Model and Set variable from combined model variable
        annInfectionsDf[c('Sim', 'GCM', 'RCP')] <- str_split_fixed(annInfectionsDf$model, ' ', 3)
        df4 <- annInfectionsDf %>% mutate(GCM = paste0(period2,"_",GCM))
        
        df <- rbind(df1, df2, df3, df4)
        tick_label <- c('GFDL-ESM2M', 'HadGEM2-ES', 'IPSL-CM5A-LR', 'MIROC5')
        # reverse the elements in levels to have the expected order of the bars in final plot
        df$GCM <- factor(df$GCM, levels = c(as.character(unique(df$GCM))))
        cols <- c("#F8766D", "#7CAE00", "#00BFC4", "#C77CFF")
        box_cols <- c("gray", cols, cols, cols)
        plot1 <- ggplot(df , aes(x = GCM, y = infections, fill = GCM)) +
                geom_boxplot(alpha=0.3, outlier.shape = NA)+
                xlab("") + ylab("Annual total infections") +
                scale_y_continuous(limits =  c(0, 16e6), breaks=c(0, 3e6, 6e6, 9e6, 12e6, 15e6))+
                scale_x_discrete(breaks=c(as.character(unique(df$GCM))),
                                 labels=c("Observed", tick_label, tick_label, tick_label))+
                geom_vline(xintercept = c(1.5,5.5,9.5), linetype = "longdash")+
                annotate("text", x = c(3.5,7.5,11.5), y = 16e6, 
                         label = c("1995_2014","2030-2049", "2080-2099"), size = 8)+
                scale_fill_manual(values=box_cols)+
                theme(text = element_text(face = "bold",size=18,colour="black"),
                      axis.text.x = element_text(angle = 45, hjust = 1),
                      panel.background = element_blank(),
                      axis.line = element_line(colour = "black"),
                      legend.position="none")
        
        rm(df, df2, df3, df4)
        
        simFolder <- file.path(projectFolder, "savedsims", paste0(year_range), RCP2)
        # load summary file with output from 1000 sims
        annInfectionsDf <- extractorRData(file.path(simFolder,"sim_summary.Rda"), 'annInfectionsDf')
        # Create Model and Set variable from combined model variable
        annInfectionsDf[c('Sim', 'GCM', 'RCP')] <- str_split_fixed(annInfectionsDf$model, ' ', 3)
        df2 <-  annInfectionsDf %>% mutate(GCM = paste0(year_range,"_",GCM))
        
        simFolder <- file.path(projectFolder, "savedsims", paste0(period1), paste0(RCP2), paste0(baseline_sero))
        # load summary file with output from 1000 sims
        annInfectionsDf <- extractorRData(file.path(simFolder,"sim_summary.Rda"), 'annInfectionsDf')
        # Create Model and Set variable from combined model variable
        annInfectionsDf[c('Sim', 'GCM', 'RCP')] <- str_split_fixed(annInfectionsDf$model, ' ', 3)
        df3 <- annInfectionsDf %>% mutate(GCM = paste0(period1,"_",GCM))
        
        simFolder <- file.path(projectFolder, "savedsims", paste0(period2), paste0(RCP2), paste0(baseline_sero))
        # load summary file with output from 1000 sims
        annInfectionsDf <- extractorRData(file.path(simFolder,"sim_summary.Rda"), 'annInfectionsDf')
        # Create Model and Set variable from combined model variable
        annInfectionsDf[c('Sim', 'GCM', 'RCP')] <- str_split_fixed(annInfectionsDf$model, ' ', 3)
        df4 <- annInfectionsDf %>% mutate(GCM = paste0(period2,"_",GCM))
        
        df <- rbind(df1, df2, df3, df4)
        # reverse the elements in levels to have the expected order of the bars in final plot
        df$GCM <- factor(df$GCM, levels = c(as.character(unique(df$GCM))))
        plot2 <- ggplot(df, aes(x = GCM, y = infections, fill = GCM)) +
                geom_boxplot(alpha=0.3, outlier.shape = NA)+
                xlab("") + 
                scale_y_continuous(limits = c(0,16e6), breaks=c(0, 3e6, 6e6, 9e6, 12e6, 15e6))+
                scale_x_discrete(breaks=c(as.character(unique(df$GCM))),
                                 labels=c("Observed", tick_label, tick_label, tick_label))+
                geom_vline(xintercept = c(1.5,5.5, 9.5), linetype = "longdash")+
                annotate("text", x = c(3.5,7.5,11.5), y = 16e6, 
                         label = c("1995_2014","2030-2049", "2080-2099"), size = 8)+
                scale_fill_manual(values=box_cols)+
                theme(text = element_text(face = "bold",size=18,colour="black"),
                      axis.text.x = element_text(angle = 45, hjust = 1),
                      panel.background = element_blank(),
                      axis.line = element_line(colour = "black"),
                      axis.title.y =element_blank(),
                      legend.position="none")
        plot <- plot_grid(plot1, plot2, labels = c("A", "B"),
                          label_size=18)
        return(plot)
}

#' plot_Seasonality() is an R function that plots the median monthly infections for a given baseline 
#' sero-prevalence, and two time periods and two Representative Concentration Pathways.
#' 
#' @param baseline_sero a string that specifies the baseline seroprevalence used in the simulation
#' @param period1 a string that specifies the first time period of the simulation
#' @param period2 a string that specifies the second time period of the simulation
#' @param RCP1 a string that specifies the first RCP scenario used in the simulation
#' @param RCP2 a string that specifies the second RCP scenario used in the simulation
#' @param year_range a string that specifies the year range of the observations. 
#' 
plot_Seasonality <- function(baseline_sero, period1 = "2030_2049", period2 = "2080_2099", 
                             RCP1 = "RCP 4.5", RCP2 = "RCP 8.5", year_range="1995_2014"){
        
        periods <- c(period1, period2)
        RCPs <- c(RCP1, RCP2)
        
        simFolder <- file.path(projectFolder, "savedsims", paste0(year_range), "Observed")
        monthInfectionsDf <- extractorRData(file.path(simFolder,"sim_summary.Rda"), 'monthInfectionsDf')
        monthInfectionsDf[c('Sim', 'GCM', 'RCP')] <- str_split_fixed(monthInfectionsDf$model, ' ', 3)
        dat1 <- monthInfectionsDf %>% 
                group_by(month, GCM) %>%
                summarise(medInf = median(infections)) 
        
        tempdat <- list()
        for (i in 1:length(periods)){
                for (j in 1:length(RCPs)){
                        
                        simFolder <- file.path(projectFolder, "savedsims", paste0(periods[i]), paste0(RCPs[j]), paste0(baseline_sero))
                        monthInfectionsDf <- extractorRData(file.path(simFolder,"sim_summary.Rda"), 'monthInfectionsDf')
                        monthInfectionsDf[c('Sim', 'GCM', 'RCP')] <- str_split_fixed(monthInfectionsDf$model, ' ', 3)
                        name <- paste(periods[i],RCPs[j])
                        tempdat[[name]] <- monthInfectionsDf %>% 
                                group_by(month, GCM) %>%
                                summarise(medInf = median(infections))
                        
                }
        }
        
        dat <- unnest(enframe(tempdat))
        
        plot <- ggplot(dat) +
                geom_line(aes(x = month, y=medInf, color = GCM))+
                geom_point(aes(x = month, y=medInf, color = GCM))+
                annotate(geom='line', x=dat1$month,y=dat1$medInf)+
                # scale_y_log10()+
                xlab("Month") + ylab("Monthly infections") +
                # scale_y_continuous(limits = c(0,35))+
                scale_x_continuous(breaks = round(seq(1, 12, by = 1), 1),
                                   labels = c("J","F","M","A","M","J","J","A","S","O","N","D"))+
                facet_wrap(~name, ncol = 2)+
                
                theme(legend.title = element_blank(), legend.position="top",
                      text = element_text(face = "bold",size=18,colour="black"),
                      panel.background = element_blank(),
                      axis.line = element_line(colour = "black"))
        return(plot)
        
}

#' seasonal_plot_Nv() is an R function that plots the median monthly mosquito numbers for a given baseline 
#' sero-prevalence, and two time periods and two Representative Concentration Pathways.

#' @param baseline_sero a string that specifies the baseline seroprevalence used in the simulation
#' @param period1 a string that specifies the first time period of the simulation
#' @param period2 a string that specifies the second time period of the simulation
#' @param RCP1 a string that specifies the first RCP scenario used in the simulation
#' @param RCP2 a string that specifies the second RCP scenario used in the simulation
#' 

seasonal_plot_Nv <- function(baseline_sero, period1 = "2030_2049", period2 = "2080_2099", 
                    RCP1 = "RCP 4.5", RCP2 = "RCP 8.5"){
        
        periods <- c(period1, period2)
        RCPs <- c(RCP1, RCP2)
        tempdat <- list()
        
        for (i in 1:length(periods)){
                for (j in 1:length(RCPs)){
                        
                        simFolder <- file.path(projectFolder, "savedsims", paste0(periods[i]), paste0(RCPs[j]), paste0(baseline_sero))
                        datDf <- extractorRData(file.path(simFolder,"sim_summary.Rda"), 'datDf')
                        # Create Model and Set variable from combined model variable
                        datDf[c('Sim', 'GCM', 'RCP')] <- str_split_fixed(datDf$model, ' ', 3)
                        name <- paste(periods[i],RCPs[j])
                        tempdat[[name]] <- datDf %>% 
                                mutate(month = month(date)) %>%
                                group_by(month, GCM) %>%
                                summarise(medNv = median(Nv)) %>%
                                # replace really small values with 1
                                mutate(medNv = ifelse(medNv<1, 1, medNv))
                                
                        rm(datDf)
                }
        }
        dat <- unnest(enframe(tempdat))
        
        plot <- ggplot(dat) +
                geom_line(aes(x = month, y=medNv, color = GCM))+
                geom_point(aes(x = month, y=medNv, color = GCM))+
                xlab("Month") + ylab("Mosquito population") +
                # scale_y_log10()+
                scale_x_continuous(breaks = round(seq(1, 12, by = 1), 1),
                                   labels = c("J","F","M","A","M","J","J","A","S","O","N","D"))+
                facet_wrap(~name, ncol = 2)+
 
                theme(legend.title = element_blank(), legend.position="top",
                      text = element_text(face = "bold",size=18,colour="black"),
                      panel.background = element_blank(),
                      axis.line = element_line(colour = "black"))
        
}

#' validation_plot_MonInf is a function that creates a validation plot showing monthly infections.
#' 
#' @param period: A string argument representing the period for which the simulation data is extracted. 
#'                This argument is required.
#' @param RCP1
#' @param RCP2
#' @param RCP3: Optional string arguments representing the three RCP scenarios 
#'              for which the simulation data is extracted. The default values 
#'              are "RCP 4.5", "RCP 8.5", and "Observed", respectively.


validation_plot_MonInf <- function(period, 
                                   RCP1 = "RCP 4.5", RCP2 = "RCP 8.5", RCP3 = "Observed"){
        
        RCPs <- c(RCP1) # removed RCP2 as Ian suggested to do so
        tempdat <- list()
        
        for (i in 1:length(RCPs)){
                simFolder <- file.path(projectFolder, "savedsims", paste0(period), paste0(RCPs[i]))
                monthInfectionsDf <- extractorRData(file.path(simFolder,"sim_summary.Rda"), 'monthInfectionsDf')
                monthInfectionsDf[c('Sim', 'GCM', 'RCP')] <- str_split_fixed(monthInfectionsDf$model, ' ', 3)
                name <- paste(period,RCPs[i])
                tempdat[[name]] <- monthInfectionsDf %>% group_by(year, month, GCM) %>%
                        summarise(medNewInf = median(infections)) %>%
                        # replace really small values with 1
                        mutate(medNewInf = ifelse(medNewInf<1, 1, medNewInf))
        }
        
        dat <- unnest(enframe(tempdat))
        
        simFolder <- file.path(projectFolder, "savedsims", paste0(period), paste0(RCP3))
        monthInfectionsDf <- extractorRData(file.path(simFolder,"sim_summary.Rda"), 'monthInfectionsDf')
        monthInfectionsDf[c('Sim', 'GCM', 'RCP')] <- str_split_fixed(monthInfectionsDf$model, ' ', 3)
        dat1 <- monthInfectionsDf %>% group_by(year, month, GCM) %>%
                summarise(medNewInf = median(infections)) %>%
                # replace really small values with 1
                mutate(medNewInf = ifelse(medNewInf<1, 1, medNewInf))
        
        yearPoints <- dat1$year+dat1$month/12
        ptsPlots <- yearPoints[seq(6, length(yearPoints), 12)]
        yearLabels <- paste0(1:20)
        
        plot1 <- ggplot(dat) +
                geom_line(aes(x = year+month/12, y=medNewInf, color = GCM))+
                xlab("Year") + ylab("Monthly infections") +
                
                facet_wrap(~name, ncol = 1)+
                geom_line(dat1, mapping = aes(x = year+month/12, y=medNewInf))+
                scale_x_continuous(breaks = ptsPlots,
                                   labels = yearLabels) +
                theme(legend.title = element_blank(), legend.position="none",
                      text = element_text(face = "bold",size=18,colour="black"),
                      panel.background = element_blank(),
                      axis.line = element_line(colour = "black"))
        
        plot2 <- ggplot(dat) +
                geom_line(aes(x = year+month/12, y=medNewInf, color = GCM))+
                xlab("Year") + ylab("Monthly infections (log)") +
                scale_y_log10()+
                facet_wrap(~name, ncol = 1)+
                geom_line(dat1, mapping = aes(x = year+month/12, y=medNewInf))+
                scale_x_continuous(breaks = ptsPlots, 
                                   labels = yearLabels) +
                theme(legend.title = element_blank(), legend.position="none",
                      text = element_text(face = "bold",size=18,colour="black"),
                      panel.background = element_blank(),
                      axis.line = element_line(colour = "black"))
        
        legend <- get_legend(plot1 + guides(color = guide_legend(nrow = 1)) + theme(legend.position = "top"))
        
        plot3 <- plot_grid(plot1, plot2, labels = c("A", "B"),
                           hjust = -4, label_size=18)
        
        plot <- plot_grid(legend, plot3, ncol = 1, rel_heights = c(.1, 1))
        
        return(plot)
}




validation_plot_sero_season <- function(period, RCP1 = "RCP 4.5", RCP2 = "RCP 8.5", RCP3 = "Observed"){
        
        RCPs <- c(RCP1) # removed RCP2 as Ian suggested to do so
        tempdat <- list()
        for (i in 1:length(RCPs)){
                
                simFolder <- file.path(projectFolder, "savedsims", paste0(period), paste0(RCPs[i]))
                datDf <- extractorRData(file.path(simFolder,"sim_summary.Rda"), 'datDf')
                # Create Model and Set variable from combined model variable
                datDf[c('Sim', 'GCM', 'RCP')] <- str_split_fixed(datDf$model, ' ', 3)
                name <- paste(period,RCPs[i])
                tempdat[[name]] <- datDf %>% group_by(days, GCM) %>%
                        summarise(medPrev = median(seroPrev))
                rm(datDf)
                
        }
        dat <- unnest(enframe(tempdat))
        
        simFolder <- file.path(projectFolder, "savedsims", paste0(period), paste0(RCP3))
        datDf <- extractorRData(file.path(simFolder,"sim_summary.Rda"), 'datDf')
        # Create Model and Set variable from combined model variable
        datDf[c('Sim', 'GCM', 'RCP')] <- str_split_fixed(datDf$model, ' ', 3)
        
        dat1 <- datDf %>% group_by(days, GCM) %>%
                summarise(medPrev = median(seroPrev))
        
        # Plot specifications
        yearPoints <- c(0.5:19.5)
        # yearPoints <- c(0.5,4.5,8.5,12.5,16.5,19.5)
        ptsPlots <- yearPoints*365
        yearLabels <- paste0(1:20)
        # yearLabels <- c("1995", "1999", "2003", "2007", "2011", "2015")

        plot1 <- ggplot(dat %>% filter(GCM != "Observed")) +
                geom_line(aes(x = days, y=medPrev, color = GCM))+
                
                xlab("Year") + ylab("Seroprevalence") +
                scale_y_continuous(limits = c(0,85))+
                scale_x_continuous(breaks = ptsPlots, 
                                   labels = yearLabels) +
                facet_wrap(~name, ncol = 1)+
                geom_line(dat1, mapping = aes(x = days, y=medPrev))+
                theme(legend.title = element_blank(), legend.position="none",
                      text = element_text(face = "bold",size=18,colour="black"),
                      panel.background = element_blank(),
                      axis.line = element_line(colour = "black"))
        
        # Proportion of infection falling in each month
        tempdat <- list()

        for (i in 1:length(RCPs)){
                        
                simFolder <- file.path(projectFolder, "savedsims", paste0(period), paste0(RCPs[i]))
                monthInfectionsDf <- extractorRData(file.path(simFolder,"sim_summary.Rda"), 'monthInfectionsDf')
                monthInfectionsDf[c('Sim', 'GCM', 'RCP')] <- str_split_fixed(monthInfectionsDf$model, ' ', 3)
                name <- paste(period,RCPs[i])
                tempdat[[name]] <- monthInfectionsDf %>% 
                                # group_by(Sim, GCM, year) %>%
                                # mutate(annInfection = sum(infections)) %>% 
                                # mutate(monthProp = 100*infections/annInfection) %>%
                                group_by(month, GCM) %>%
                                summarise(avgmonthProp = mean(infections))
                        

        }
        
        dat <- unnest(enframe(tempdat))
        
        simFolder <- file.path(projectFolder, "savedsims", paste0(period), paste0(RCP3))
        monthInfectionsDf <- extractorRData(file.path(simFolder,"sim_summary.Rda"), 'monthInfectionsDf')
        monthInfectionsDf[c('Sim', 'GCM', 'RCP')] <- str_split_fixed(monthInfectionsDf$model, ' ', 3)
        dat2 <- monthInfectionsDf %>% 
                # group_by(Sim, GCM, year) %>%
                # mutate(annInfection = sum(infections)) %>% 
                # mutate(monthProp = 100*infections/annInfection) %>%
                group_by(month, GCM) %>%
                summarise(avgmonthProp = median(infections))
        
        plot2 <- ggplot(dat) +
                geom_line(aes(x = month, y=avgmonthProp, color = GCM))+
                geom_point(aes(x = month, y=avgmonthProp, color = GCM))+
                # scale_y_log10()+
                xlab("Month") + ylab("Monthly Infections") +
                # scale_y_continuous(limits = c(0,35))+
                scale_x_continuous(breaks = round(seq(1, 12, by = 1), 1),
                                   labels = c("J","F","M","A","M","J","J","A","S","O","N","D"))+
                
                facet_wrap(~name, ncol = 1)+
                annotate(geom='line', x=dat2$month,y=dat2$avgmonthProp)+
                theme(legend.position="none",
                      text = element_text(face = "bold",size=18,colour="black"),
                      panel.background = element_blank(),
                      axis.line = element_line(colour = "black"))
        
        legend <- get_legend(plot1 + guides(color = guide_legend(nrow = 1)) + theme(legend.position = "top"))
        
        plot3 <- plot_grid(plot1, plot2, labels = c("A", "B"),
                           hjust = -3, label_size=18)
        
        plot <- plot_grid(legend, plot3, ncol = 1, rel_heights = c(.1, 1))
        
        return(plot)
        
        
}
