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
smc_plot_R = function(mat,days,quant=c(0.95,0.90,0.75,0.5),col = "#08769c",
                      title="",trunc=NULL,mat2=NULL,mat3=NULL,pnts=NULL,pnts2=NULL,pnts3=NULL,
                      lin=NULL,lin3=NULL,addone=TRUE,
                      nam0="R",nam1="Ins",nam2="Hosp(B)",nam3="Test(R)",nam4="Seed(G)",sc=0.025,sc2=0.0025,sc3=0.001,
                      end_date=end_date,test_date=NULL,pred_date=NULL,log="")
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
  show(dim(mat_melt))
  show(dim(lin3))
  q = ggplot(data=mat_melt,aes(x = dates, y = value)) + ggtitle(title) + 
    stat_lineribbon(aes(y = value), .width = quant, color = "black") +  
    scale_fill_brewer(palette="Blues") + 
    scale_x_date(date_breaks = "2 month", date_labels = "%b") +
    scale_y_continuous(name=nam0) +
    scale_x_date(date_breaks = "2 month", date_labels = "%b")+
    theme(axis.text=element_text(size=20),axis.title=element_text(size=22,face="bold"),
        text=element_text(size=18),legend.text = element_text(size=18))
  if(addone)
    q = q + geom_hline(yintercept=1,size=1.5, linetype='dashed', col = 'black')
  if(!is.null(test_date))
    q = q + geom_vline(xintercept=test_date,size=1.5, linetype='dashed', col = 'green')
  if(!is.null(pred_date))
    q = q + geom_vline(xintercept=pred_date,size=1.5, linetype='dashed', col = 'green')
  if(!is.null(mat2))
  {
    quant2 = sort(unique(c((1-quant)/2,1-(1-quant)/2)))
    show(quant2)
    ci = apply(mat2,1,quantile,probs=quant2,na.rm=TRUE)
    if(!is.null(trunc))
      ci = pmin(ci,trunc)

    d = data.frame(days=days,ci=t(ci))
    names(d)= c("days",paste("ci",1:nrow(ci),sep=""))
    q = q + geom_path(aes(x = days, y = ci3),data=d,color="red")
    q = q + geom_path(aes(x = days, y = ci6),data=d,color="red")
    q = q + geom_abline(intercept=1.00,slope=0.0,color="red")
  }
  if(!is.null(mat3))
  {
    ci = mat3
    if(!is.null(trunc))
      ci = pmin(ci,trunc)
    
    d = data.frame(days=days,ci=(ci))
    names(d)= c("days",paste("ci",1:ncol(ci),sep=""))
    q = q + geom_path(aes(x = days, y = ci1),data=d,color="red",linetype=2)
    q = q + geom_path(aes(x = days, y = ci2),data=d,color="red",linetype=1)
    q = q + geom_path(aes(x = days, y = ci3),data=d,color="red",linetype=2)
  }
  if(!is.null(lin))
  {
    d = data.frame(days=days,Rtrue=lin)
    q = q + geom_path(aes(x=days,y=lin),data=d,color="red",linetype=1)
  }
  if(!is.null(lin3))
  {
   q = q + scale_fill_manual(values = c("skyblue3", "black")) 
    lin3 = pmin(lin3,trunc)
   for(i in 1:3)
   {
    d = data.frame(days=days,Rtrue=lin3[,i])
    q = q + geom_path(aes(x=days,y=Rtrue),data=d,color="red",linetype=min(i,2))
   }
  }
  if(!is.null(pnts))
  {
    q = q + geom_point(data=pnts,mapping=aes(date,N*sc),color="blue",size=0.2)
    q = q + scale_y_continuous(sec.axis = sec_axis(~./sc, name = nam2))
  }
  if(!is.null(pnts2))
  {
    q = q + geom_point(data=pnts2,mapping=aes(date,TP*sc2),color="red",size=0.2)
    q = q + scale_y_continuous(sec.axis = sec_axis(~./sc2, name = paste(nam2,nam3)))
  }
  if(!is.null(pnts3))
  {
    q = q + geom_point(data=pnts3,mapping=aes(date,Seed*sc3),color="green",size=0.2)
    q = q + scale_y_continuous(sec.axis = sec_axis(~./sc3, name = paste(nam2,nam3,nam4)))
  }
  q
}
