#' Plotting routine for estimates of R. Make plots for the given set of quantiles
#'
#' @param mat Matrix of R-simulations, rows corresponding to days, columns to particles
#' @param days Vector of dates, length should correspond to number of rows in 'mat'
#' @param quant Quantiles to be plotted, default is c(0.95,0.9,0.5,0.01)
#' @param col Colour code used for the different quantiles. Default is col = "#08769c" giving different blue colors
#' @param title Optional title to include to the plot. Default is NULL
#' @param trunc If included, plots are truncated from above along the y-axis for this value
#' @param mat2 An additional alternative set of R-simulations. If included, 50 percent quantiles are added (with green color) on the plot
#' @param end_date Last date to be plotted.
#' @param test_date Time for which test data are included. If given, a vertical red line is given at this date
#' @param pred_date Time for which prediction starts. If given, a vertical red line is given at this date
#' @return An object of class gg / ggplot which can be plotted
#' @export
smc_plot_data = function(mat,days,datD,quant=c(0.95,0.90,0.75,0.5),col = "#08769c",
                      title="",trunc=NULL,nam="Hospital", 
                      start_date=start_date,end_date=end_date,test_date=NULL,pred_date=NULL,log=FALSE)
{
  require(ggplot2)
  require(tidybayes)

  if(!is.null(trunc))
    mat = pmin(mat,trunc)
  #Add dates
  mat<-data.table(mat)
  mat[,dates:=days]
 #Wide to long format
  mat_melt<-melt(mat,id.vars="dates")
  mat_melt<-mat_melt[dates<=end_date,]
  q = ggplot(data=mat_melt,aes(x = dates, y = value)) + ggtitle(title) + 
    stat_lineribbon(aes(y = value), .width = quant, color = "black") +  
    scale_fill_brewer(palette="Blues") +
    scale_x_date(date_breaks = "2 month", date_labels = "%b") +
    theme(text=element_text(size=18))+
    scale_x_date(date_breaks = "2 month", date_labels = "%b") +
    theme(axis.text=element_text(size=20),axis.title=element_text(size=22,face="bold"),
        text=element_text(size=18),legend.text = element_text(size=18))
 if(log)
    q = q + scale_y_continuous(name=nam,trans='log10')
  else
    q = q + scale_y_continuous(name=nam)
  if(nam=="Hospital")
    q = q + geom_point(datD,mapping=aes(date,N),color="red",size=0.5)
  else if(nam=="Test")
    q = q + geom_point(datD,mapping=aes(date,Positive),color="red",size=0.5)
  
  if(!is.null(test_date))
    q = q + geom_vline(xintercept=test_date,size=1.5, linetype='dashed', col = 'green')
  if(!is.null(pred_date))
    q = q + geom_vline(xintercept=pred_date,size=1.5, linetype='dashed', col = 'green')
 q
}
