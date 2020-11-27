
bilinear_model = function(X, D = diag(ncol(X)), lambda, ncomp, verbose = F, init = c("SVD", "random")){
  if (nrow(D) == ncol(X) && nrow(D) == nrow(X)) stop("D has the wrong dimensions")
  p = ncol(X)
  n = nrow(X)
  init = match.arg(init)
  if (init == "random"){
    A = matrix(runif(n*ncomp), nrow = n, ncol = ncomp)
    B = matrix(runif(p*ncomp), nrow = p, ncol = ncomp )
  }
  if (init == "SVD" ){
    s = svd(X, nu = ncomp, nv = ncomp)
    A = s$u
    B = s$v
  }

  for (i in 1:ncomp){
    if (verbose == TRUE) print(paste("Compoenent: ", i, "Lambda: ", lambda ))
    a = A[, i, drop = F]
    X_tilde = X - A[, -i, drop = F] %*%t(B[, -i, drop = F])  #4
    u = Variable(p)
    objective = Minimize( 0.5 * cvxr_norm( t(X_tilde) %*% a - t(D) %*% u ))
    constraint = list(cvxr_norm(u, p = "inf") <= lambda)

    problem = Problem(objective,constraint)
    sol = solve(problem,)#5
    u_hat = sol$getValue(u)
    B[,i] = t(X_tilde)%*%a - t(D) %*% u_hat
    A[, i] = X_tilde%*%B[,i]/norm(X_tilde%*%B[,i], type = "2")
  }
  rss = sum((X - A%*%t(B))^2 )
  obj = list(A = A, B = B, rss = rss ,
             pearson = cor(as.numeric(X), as.numeric(A%*%t(B))),
             evar = 1 - (rss/sum(X^2)) )
  attr(obj, "class") = "factor_model"
  return(obj)
}
