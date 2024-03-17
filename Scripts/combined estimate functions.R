#some functions for later

Ne_to_r <- function(Ne, S){
  r <- (-69*S^2 + sqrt(10000 * (S^4) * (Ne^2) + 4761 * (S^4) - 248400*(S^3)*Ne) + 1800*S*(Ne^2))/( 1800*(S^2)*(Ne^2))  
  return(r)
}

Ne_to_r_star <- function(Ne, S){
  r <- Ne_to_r(Ne, S)
  r <- r - 1/S
  return(r)
}