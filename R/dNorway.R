#' Example data set for Covid-19 data in Norway
#'
#' A dataset containing the necessary data to run the SMC algorithm. The hospital data, the number of positives in the test data and the initial population sizes are real data, others are (due to confidentiality) simulated so that they resmble the real data. The variables are as follows:
#'
#' \itemize{
#'   \item hosp: A data.table with two columns, 'date' and 'N' where 'N' is the number of new hospital admissions at day 'date'
#'   \item test: A data.table with three columns, 'date', 'Tot', 'Positive' with 'Tot' being total tests and 'Positive' being number of positives
#'   \item seed: A datatable with 3 columns, 'location_code' (location of region), 'day' 'n' (number of imported cases). Note that 'day' here is a number indicated the number of days from start_date!
#'   \item mob_mat: A list with each element containing mobility data for 6 hour intervals, so mob_mat[[(i-1)*4+j]] correspond to day i, period j on that day. 
#'   Each element in the list is a data.table with three columns, the two first giving regions from and to the last giving numbers of individuals that are moved.
#'   \item pop: A data table of dimension nreg x 7 giving the initial compartmental values for each region
#'   \item start_date: The starting data for the 'seed' and 'mob_mat' data. Should be in correspondence with the starting data used when running the SMC algorithm
#' }
#'
#' @docType data
#' @keywords datasets
#' @name dNorway
#' @usage data(dNorway)
#' @format A list of 5 objects, each explained above though
NULL
