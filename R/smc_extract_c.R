#' Extract c_symb_incidence+c_asymp_incidence summed by day from a (list of) data tables
#'
#' @param 'particle' A data table containing at least c_symb_incidence and c_asymp_incidence
#' @param 'particles' A list of particles, each a data table containing at least c_symb_incidence and c_asymp_incidence
#' @return A vector/matrix with rows corresponding to days and if calling smc_extract_list, a column corresponding to particles
smc_extract_c = function(particle)
{
  particle[,.(cc=sum(c_symp_incidence+c_asymp_incidence)),by=day]$cc
}

smc_extract_c_list = function(particles)
{
  B = length(particles)
  matrix(unlist(lapply(particles,smc_extract_c)),ncol=B)
}