#' @import data.table 
#' @export
#'
get_seeding = function(dir="/mn/sarpanitu/ansatte-u2/geirs/prj/Covid-19/data_reg/",
                       file = "seeding.RDS",
                       first_seeding_date= "2020-02-17"){
        #Creating the seeding matrix
 inf.data <- readRDS(paste(dir,file,sep=""))
 setDT(inf.data)
 infe.seed <- inf.data
 infe.seed <- infe.seed[date >first_seeding_date,]
 infe.seed[, day:=round(lubridate::interval(first_seeding_date, date)/lubridate::days(1))]
 infe.seed <- infe.seed[,c("municip_code", "day", "N")]
 colnames(infe.seed) <- c( "location_code", "day", "n")
 infe.seed <- infe.seed[ infe.seed$n!= 0, ]
 #  infe.seed <- infe.seed[, .(n=sum(n)), by=.(location_code, day)]
 return(infe.seed)
}

