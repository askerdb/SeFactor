simulate_factor_model = function(n_art = 40, p_art = 50,k_art = 2, sigma = 1){
  A_art = mvrnorm(n_art, mu = rep(0, k_art,), diag(k_art))
  B_art = sapply(1:k_art, function(x) rbinom(p_art, 1,  prob = 0.5))
  E_art = mvrnorm(n = n_art, mu = rep(0, p_art), Sigma = sigma^2 * diag(p_art))
  return(list(X = A_art%*%t(B_art), A_art = A_art, B_art = B_art))
}
