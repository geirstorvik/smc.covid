#' Plotting routine for estimates of R. Make plots for the given set of quantiles
#'
#' @param 'mat' Matrix of R-simulations, rows corresponding to days, columns to particles
#' @param 'days' Vector of dates, length should correspond to number of rows in 'mat'
#' @param 'quant' Quantiles to be plotted, default is c(0.95,0.9,0.5,0.01)
#' @param 'col' Colour code used for the different quantiles. Default is col = "#08769c" giving different blue colors
#' @param 'title' Optional title to include to the plot. Default is NULL
#' @param 'trunc' If included, plots are truncated from above along the y-axis for this value
#' @param 'mat2' An additional alternative set of R-simulations. If included, 50 percent quantiles are added (with green color) on the plot
#' @param 'end_date' Last date to be plotted.
#' @param 'test_date' Time for which test data are included. If given, a vertical red line is given at this date
#' @param 'pred_date' Time for which prediction starts. If given, a vertical red line is given at this date
#' @return An object of class gg / ggplot which can be plotted
smc_plot_cmp_median_R = function(mat,days,col = "#08769c",
                      title="",
                      end_date=end_date,test_date=NULL,pred_date=NULL)
{
  require(ggplot2)
  require(tidybayes)

  col = c("black","red","blue","green")
  
  #Add dates
 mat<-data.table(mat)
  mat[,dates:=days]
  #Wide to long format
  mat_melt<-melt(mat,id.vars="dates")
  show(dim(mat_melt))
  mat_melt<-mat_melt[dates<=end_date,]
  show(dim(mat_melt))
  q = ggplot(mat_melt, aes(x = dates, y = value)) + 
    geom_line(aes(color = variable),size=2) + scale_color_manual(values=c("black", "red","blue","green"))  + 
    scale_x_date(date_breaks = "1 month", date_labels = "%b") + 
      theme(axis.text=element_text(size=14),axis.title=element_text(size=14,face="bold"),legend.text=element_text(size=14))
  q
}

