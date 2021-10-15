#Script for running the SMC algorithm on a dataset based on real COVID-19 data from Norway
#New hospital incidences and daily reported positive PCR tests are real data
#Imported cases (seeds) are simulated from a Poisson distribution with intensities fitted from real data
#Total number of tests are simulated from a Poisson distribution with intensities runnin means of the real data.
#Mobility matrices are assumed equal for all days but different for the 4 6 hours intervals, mean over the whole period
library(rstudioapi)
setwd(dirname(getActiveDocumentContext()$path)) 
library(smc.covid)
#Model spesifications
dyn_model = "AR_mu_0"               #AR model with mu=0
#dyn_model = "AR"
#dyn_model = "RW"
#dyn_model = "CP"

#Number of particles, 200 here for relatively fast computation and testing of whether it works
#Takes about ?? on a linux server using 4 cores
#Increasing to 2000 should give reasonable estimates of the R process for this data set
#In the paper we have used 20000 which takes 4-5 hours on a linux server with 50 cores
B = 2000
N_cores = 25

#Directory for data and storing results
mdir = "./"
resdir = paste(mdir,"test_run_",dyn_model,"_B_",B,"/",sep="")

#Run settings
RUN = FALSE               #Will run the SMC algorithm and store the results in subdirectory of resdir
PLOT = TRUE              #With RUN=FALSE, plotting can be performed on previously stored results

#Seed
inseed = 44

#DATES
start_date<-as.Date("2020-02-17")
end_date<-as.Date("2021-07-01")                 #End date of observations
end_date2 = end_date+lubridate::days(21)        #End date for prediction
ndays = end_date-start_date                     #Number of days with data
ndays2 = end_date2-start_date                   #Total number of days
days = start_date+lubridate::days(0:ndays2)
npred = as.numeric(end_date2-end_date)         #Number of days for forecasting

#Read data
data(dNorway)

#Remove data after end_date
hospD = dNorway$hosp
hospD$N[hospD$date>end_date] = NA
testD = dNorway$test
testD$Positive[testD$date>end_date] = NA


#For prediction, putting seeds (imported cases), and mobility equal to values three weeks 
loc.seed = dNorway$seed
loc.seed = loc.seed[day<=ndays]
x = loc.seed[day>(ndays-21)]
x$day = x$day+21
for(day in 1:npred)
{
  x = loc.seed[day==(ndays+day-21)]
  
  loc.seed[day==(ndays+day)] = loc.seed[day==(ndays+day-21)]
}
loc.seed = rbind(loc.seed,x)
n = 4*ndays
mob_mat = dNorway$mob_mat[1:n]
for(i in 1:21)
  for(j in 1:4)
  {
    mob_mat[n+4*(i-1)+j] = mob_mat[n+4*(i-1-21)+j]
  }
#For prediction, total tests are averages of values in previous weeks
testD2 = testD
for(d in 1:21)
{
  testD2$Tot[335+d] = as.integer(mean(testD2$Tot[335+d-c(1:7)]))
}

if(RUN)    ## Run SMC algorithm
{
  dir.create(resdir)
  res.smc<-smc_loglik_mem(mu_R0 = 3.3,sd_R0 = 0.1,mu_AMP = 2.8,sd_AMP = 0.0,
                    param_dynamicR = list(alpha=2.4,beta=0.28,phi.a=0.5,phi.b=9,a0=0.5,mu0=1,kappa0=16),
                    dyn_model=dyn_model,
                    lambda = c(8.87,7.60),lambda_size=c(3.40,3.34),lambda_ch=as.Date("2020-08-01"), 
                    prob_h = c(0.0378,0.0230,0.0209,0.0180,0.0146,0.0185,0.0183,0.0208,0.0232,0.0221,0.0230,0.0224,0.0231, 0.0141,0.0092,0.0098),
                    q_test=c(0,0.05,0.080,0.16,0.4,0.3,0.01),
                    change_prob=c(as.Date("2020-05-01"),as.Date("2020-06-01"),as.Date("2020-07-01"),
                                  as.Date("2020-08-01"),as.Date("2020-09-01"),as.Date("2020-10-01"),as.Date("2020-11-01"),
                                  as.Date("2020-12-01"),as.Date("2021-01-01"),as.Date("2021-02-01"),as.Date("2021-03-01"),as.Date("2021-04-01"), as.Date("2021-05-01"),as.Date("2021-06-01"),
                                  as.Date("2021-07-01")),
                    alpha_hosp =8,alpha_test=8, pi_test = c(-1.0173862969, 0.0001318322),
                    start_date = start_date,dynR_date = as.Date("2020-03-08"),end_date = end_date2,
                    hosp_data = hospD,test_data=testD2,   
                    loc.seed=loc.seed,mob_mat=mob_mat,pop=dNorway$pop,
                    resdir=resdir,
                    B = B,unit.samp = 3,resamp.frac = 1.1,smo.lag = 25,
                    parallel = TRUE,sim.data=TRUE,
                    VERBOSE=TRUE,output_latent = FALSE,seed = inseed,N_cores=N_cores)
  #Store output in resdir
  #NOTE: The filtered estimates are stored dynamically in separate files in the resdir directory, see below
  saveRDS(res.smc,paste(resdir,"Norway_res.RDS",sep=""))
}
if(PLOT)
{
  if(!exists("res.smc$"))
    res.smc = readRDS(paste(resdir,"Norway_res.RDS",sep=""))
  
  #Plot fixed lag predictions of the R-process
  matR = as.matrix(read.table(paste(resdir,"smc_dynamic_R_lag.txt",sep="")))
  pdf(paste(resdir,"Res_R_nonsmooth_Norway_",dyn_model,"_B",B,".pdf",sep=""),height=5,width=15)
  print(smc_plot_R(matR,days,end_date=end_date2,test_date=as.Date("2020-08-01"),pred_date=end_date,trunc=3.1))
  dev.off()
  #Plot fixed lag predictions of the R-process (weekly smoothed on log-scale)
  matR.s = smc_smooth_week_log(matR)
  pdf(paste(resdir,"Res_R_Norway_",dyn_model,"_B",B,".pdf",sep=""),height=5,width=15)
  print(smc_plot_R(matR,days,end_date=end_date2,test_date=as.Date("2020-08-01"),pred_date=end_date,trunc=3.1))
  dev.off()
  
  
  # Plot new incidences
  matC.s = smc_smooth_week(as.matrix(read.table(paste(resdir,"smc_dynamic_C_lag.txt",sep=""))))
  matC.s[1:10,] = NA
  pdf(paste(resdir,"Res_C_Norway_",dyn_model,"_B",B,".pdf",sep=""),height=5,width=15)
  print(smc_plot_C(matC.s,days,end_date=end_date2,test_date=as.Date("2020-08-01"),pred_date=end_date,trunc=2000))
  dev.off()
  
  # Plot predicted hospital data, including actual observations
  n = hospD$date[1]-start_date
  hospD2 = dNorway$hosp
  hospD2$N[hospD$date>end_date] = NA
  #hospD2 = hospD
  for(i in (n-1):0)
    hospD2 = rbind(list(start_date+i,NA),hospD2) 
  pdf(paste(resdir,"Res_Hosp_Norway_",dyn_model,"_B",B,".pdf",sep=""),height=5,width=15)
  print(smc_plot_data(res.smc$hospD.sim,days,end_date=end_date2,test_date=as.Date("2020-08-01"),pred_date=end_date,
                      datD=hospD2[hospD2$date<=end_date],trunc=60))
  dev.off()
  
  # Plot predicted test data, including actual observations
  n = testD$date[1]-start_date
  testD2 = testD
  for(i in (n-1):0)
    testD2 = rbind(list(start_date+i,NA,NA),testD2) 
  #Fill in observations in forcasting interval
  testD2$Positive[502:522] = dNorway$test$Positive[336:356]
  pdf(paste(resdir,"Res_Test_Norway_",dyn_model,"_B",B,".pdf",sep=""),height=5,width=15)
  print(smc_plot_data(res.smc$testD.sim,days,start_date=start_date,end_date=end_date2,test_date=as.Date("2020-08-01"),pred_date=end_date,
                      nam="Test",datD=testD2[testD2$date<=end_date],trunc=1500))
  dev.off()
 
  #Plotting dynamic parameter estimates, which parameters depend on model
  #Also print out final estimates (and credibility intervals)
  npar = 2
  nam = c("a","sigma")
  if(dyn_model=="RW")
  {
    npar = 1
    nam = "sigma"
  }
  if(dyn_model=="AR")
  {
    npar = 3
    nam = c("a","sigma","mu")
  }
  if(dyn_model=="AR_mu_0")
  {
    npar = 2
    nam = c("a","sigma")
  }
  if(dyn_model=="CP")
  {
    npar = 2
    nam = c("psi","sigma")
  }
  for(j in 1:npar)
  {
    matp = as.matrix(read.table(paste(resdir,"smc_SdR.txt",j,sep="")))
    pdf(paste(resdir,"Res_par_Norway_",dyn_model,"_B",B,"_par",j,".pdf",sep=""),height=5,width=15)
    print(smc_plot_data(matp[1:ndays,],days[1:ndays],end_date=end_date2,test_date=as.Date("2020-08-01"),trunc=2000,nam=nam[j]))
    dev.off()
    show(nam[j])
    show(round(c(mean(matp[ndays,]),quantile(matp[ndays,],probs=c(0.025,0.975))),3))
  }
}

