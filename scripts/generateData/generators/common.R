library(spatstat)

get_kde_sample <- function(
    init_values = ? numeric, sample_size = ? integer
){
  kde_den <- stats::density(init_values, n=5000)
  tt <- quantile(kde_den, probs = stats::runif(sample_size))
  return(as.numeric(tt))
}