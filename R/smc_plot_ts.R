#' @export
smc_plot_ts = function(mat,days,title="")
{
 require(ggplot2)
 require(tidybayes)

 #Add dates
 mat<-data.table(mat)
 mat[,dates:=days]
 #Wide to long format
 mat_melt<-melt(mat,id.vars="dates")
 q <- ggplot(data=mat_melt,aes(x = dates, y = value,colour=variable))+ 
  geom_line() + ggtitle(title)
#  geom_path(data=mat_melt,aes(x = dates, y = value),colour=variable)+
#  scale_color_manual(values = c("black", "red")) 
 q
}
