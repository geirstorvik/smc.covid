#' New run_model for smc
#'
#' @param seeding -> Data frame including importations. 3 columns: location_code, day (day 1 correspond to start date), n
#' @param se1e2iiar_pop-> Data frame with the locations. Columns: location_code,S,E1,E2,I,Ia,R.  It defines the number of locations
#' @param N_sim->Number of simulations per run (spread package returns the average of all of them)

smc_run_model <- function(mobility_matrices,
                      start_date,
                      end_date,
                      R_0=2,
                      Reffs=c(), # One national number, or a matrix which is number of locations times number of change points
                      AMP_factor=1.5,
                      reff_change_dates=c(),
                      latent_period=2.0,
                      presymptomatic_period = 2.0,
                      infectious_period = 3.0,
                      presymptomatic_relative_infectiousness = 1.30,
                      asymptomatic_prob = 0.4,
                      asymptomatic_relative_infectiousness = 0.1,
                      N_sim=100,
                      seeding=NULL,
                      se1e2iiar_pop=NULL,
                      parallel=T,
		      inputSeed = NULL){

  if(is.null(se1e2iiar_pop)){
    stop("No se1e2iiar_pop file found")
  }
  
  
  
  if (is.null(seeding)){
    print("Note: No seeding file provided.")
  }  
   
  
  N_days <- lubridate::interval(start_date, end_date) / lubridate::days(1) + 1

  mobility_matrices <- mobility_matrices[1:(N_days*4)]
  mun <- se1e2iiar_pop$location_code
  Rs <- NULL
  loc <- NULL
  days <- NULL
  times <- NULL

    if(!is.null(Reffs)){
     if(end_date >= reff_change_dates[1]){
       end_r0 <- reff_change_dates[1]
     }
    }  
    else{
      end_r0 <- end_date + 1
    }
    n1 <- lubridate::interval(start_date, end_r0) / lubridate::days(1)
    Rs <- c(Rs, rep(rep(R_0, 4 * n1), length(mun)))
    loc <- c(loc, rep(mun, each = 4 * n1))
    days <- c(days, rep(rep(1:n1, each = 4), length(mun)))
    times <- c(times, rep(c(0, 6, 12, 18), n1 * length(mun)))
    if(!is.null(Reffs)){
     if(end_date >= reff_change_dates[1]){
       for(i in 1:length(Reffs)){
          Reff <- Reffs[i]
          if(i == length(Reffs)){
            end_reff_date <- end_date + 1
          } else{
            end_reff_date <- reff_change_dates[i+1]
          }
          n2 <- lubridate::interval(reff_change_dates[i], end_reff_date) / lubridate::days(1)
          Rs <- c(Rs, rep(rep(Reffs[i], 4* n2), length(mun)))
          loc <- c(loc, rep(mun, each = 4 * n2))
          days <- c(days, rep(rep((n1 + 1):(n2 + n1), each = 4), length(mun)))
          times <- c(times, rep(c(0, 6, 12, 18), n2 * length(mun)))
          n1 <- n1 + n2
        }
      }
    }
  

  betas <-spread:::se1e2iiaR_calculate_beta_from_r0(
                     r0 = Rs,
                     a2 = 1/presymptomatic_period,
                     gamma = 1/infectious_period,
                     presymptomaticRelativeInfectiousness = presymptomatic_relative_infectiousness,
                     asymptomaticProb =  asymptomatic_prob,
                     asymptomaticRelativeInfectiousness = asymptomatic_relative_infectiousness)
		     
  betas <- data.table("location_code" = loc, "day" = days, "time" = times, "beta" = betas)
  betas <- betas[betas$day <= N_days, ]
  results <- list()
  seeding_period <- lubridate::interval(start_date, end_date) / lubridate::days(1)

  inner_run_model_split <-function(i){
    if(is.null(inputSeed)){
      inputSeed <- i + as.numeric(Sys.time())
    }
    set.seed(inputSeed)
    new_seeding <- copy(seeding)
    new_seeding$n <- new_seeding$n + rpois(length(new_seeding$n), lambda=new_seeding$n * (AMP_factor-1))
    start_se1e2iiar_pop <- copy(se1e2iiar_pop)
    res <- spread:::asymmetric_mobility_se1e2iiar(
      se1e2iiar_pop = start_se1e2iiar_pop,
      mobility_matrix = mobility_matrices[1:(seeding_period*4)],
      dynamic_seeds =  new_seeding, 
      betas = betas[betas$day <= seeding_period, ],
      latent_period = latent_period,
      infectious_period = infectious_period,
      presymptomatic_period = presymptomatic_period,
      presymptomatic_relative_infectiousness =  presymptomatic_relative_infectiousness,
      asymptomatic_prob = asymptomatic_prob,
      asymptomatic_relative_infectiousness = asymptomatic_relative_infectiousness,
      N = 1,
      inputSeed = inputSeed)
    
    res[, date:=start_date+day]
    set.seed(Sys.time())
    return(copy(res))
  }
  
  if(parallel){
    results <- parallel::mclapply(1:N_sim, inner_run_model_split, mc.cores=parallel::detectCores())
  }else{
    results <- lapply(1:N_sim, inner_run_model_split)
  }
  
  return(results)
  
}
