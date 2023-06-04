#' A function to specify ggplot options
#'
#' This function stores all my favoured plotting options and themes for
#' ggplot plots.
#'
#' @param opt Name of option to use in plot.
#' @return ggplot theme options object.
#'
#' @author Richard T. Gray, \email{Rgray@kirby.unsw.edu.au}
#'
#' @examples
#' PlotOptions()
#'
#' @export
#' @import ggplot2
#'
PlotOptions <- function(opt=c("main")) {
  
  
  opt <- match.arg(opt)
  
  # Baseline theme for plot variables
  main <- theme_bw() + theme(text = element_text(face = "bold",size=12,colour="black"),
                             axis.text.x = element_text(face = "plain",size=10,colour="black"),
                             axis.text.y = element_text(face = "plain",size=10,colour="black"),
                             axis.line.x = element_line(colour = "black"),
                             axis.line.y = element_line(colour = "black"),
                             axis.ticks = element_line(colour="black"),
                             legend.position = "top",
                             legend.background = element_rect(),
                             legend.key = element_blank(),
                             panel.grid.major = element_blank(),
                             panel.grid.minor = element_blank(),
                             panel.background = element_blank(),
                             panel.border = element_rect(colour = "black"),
                             axis.line = element_line(colour = "black"),
                             plot.title=element_text(size=12, face="bold"),
                             strip.background = element_blank()
  )
  
  # Return theme
  plotOpts <- switch(match.arg(opt),
                     main=main)
  
  return(plotOpts)
  
}