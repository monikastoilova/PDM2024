# Libraries
library(evd)
library(mev)
library(MASS)
library(doParallel)
library(foreach)

# ==============================================================================
# Section "Hybrid distributions"
# 1. Hybrid uniform-Weibull with parameters a, b, sh, sc, mix_par:
#    U(a, b) on [a, b] and Weibull(sh, sc) on (b, infty) 
#    with mixing parameter mix_par
# 2. Hybrid uniform-GP with parameters a, b, sh, sc, mix_par:
#    U(a, b) on [a, b] and GP(sc, sh) on (b, infty) 
#    with mixing parameter mix_par


# ------------------------------------------------------------------------------
r_hybrid_unif_weibull = function(n, a, b, sh, sc, mix_par){
  # function generating n data points from a hybrid 
  # uniform-Weibull distribution
  
  x = numeric(n)
  for (i in 1:n) {
    if (runif(1) <= mix_par) {
      x[i] = runif(1, min = a, max = b)
    }else{
      x[i] = b + rweibull(1, shape = sh, scale = sc)
    }
  }
  return(x)
}
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
r_hybrid_unif_gp = function(n, a, b, sh, sc, mix_par){
  # function generating n data points from a hybrid 
  # uniform-GP distribution

  x = numeric(n)
  for (i in 1:n) {
    if (runif(1) <= mix_par) {
      x[i] = runif(1, min = a, max = b)
    } else {
      x[i] = evd::rgpd(1, loc = b, scale = sc, shape = sh)
    }
  }
  return(x)
}
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
d_hybrid_unif_weibull_scalar = function(x, a, b, sh, sc, mix_par){
  # function computing the density of a hybrid uniform-Weibull 
  # distribution at x (scalar)
  if(x < a){
    return(0)
  }else if(x >= a && x <= b){
    return(mix_par * dunif(x, a, b))
  }else if(x > b){
    return((1 - mix_par) * dweibull(x - b, shape = sh, scale = sc))
  }
}
d_hybrid_unif_weibull = Vectorize(d_hybrid_unif_weibull_scalar,"x")
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
d_hybrid_unif_gp_scalar = function(x, a, b, sh, sc, mix_par){
  # function computing the density of a hybrid uniform-GP 
  # distribution at x (scalar)
  if(x < a){
    return(0)
  }else if(x >= a && x <= b){
    return(mix_par * dunif(x, 0, 1))
  }else if(x > b){
    return((1 - mix_par) * mev::dgp(x, loc = b,shape = sh, scale = sc))
  }
}
d_hybrid_unif_gp = Vectorize(d_hybrid_unif_gp_scalar,"x")
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
neg_log_likelihood_hybrid_unif_weibull = function(params, x){
  # function computing the negative log-likelihood given a 
  # vector of observations x from a hybrid uniform-Weibull 
  # distribution with parameter vector params = (a, b, sh, sc, mix_par)
  
  return(-sum(log(d_hybrid_unif_weibull(x,params[1],params[2],params[3],params[4],params[5]))))
}
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
neg_log_likelihood_hybrid_unif_gp = function(params, x){
  # function computing the negative log-likelihood given a 
  # vector of observations x from a hybrid uniform-GP 
  # distribution with parameter vector params = (a, b, sh, sc, mix_par)
  
  return(-sum(log(d_hybrid_unif_gp(x, params[1], params[2], params[3], 
                                   params[4], params[5]))))
}
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
fit_hybrid_unif_weibull = function(x, params_init = c(0, 1, 0.5, 1, 0.5)){
  # function that fits a hybrid uniform-Weibull model 
  
  fscale = neg_log_likelihood_hybrid_unif_weibull(params_init, x)
  opt_result = optim(par = params_init, 
                     neg_log_likelihood_hybrid_unif_weibull,
                     x = x, 
                     method = "Nelder-Mead", 
                     control = list(maxit = 1000, fnscale = fscale)
                     )
  res = list()
  res$estimates = opt_result$par
  names(res$estimates) = c("U:lower bound","U:upper bound",
                           "W:shape","W:scale",
                           "mixing parameter")
  res$neg_LL = opt_result$value
  return(res)
}
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
fit_hybrid_unif_gp = function(x, params_init = c(0, 1, 0.5, 1, 0.5)){
  # function that fits a hybrid uniform-Weibull model 
  
  fscale = neg_log_likelihood_hybrid_unif_gp(params_init, x)
  opt_result = optim(par = params_init, 
                     neg_log_likelihood_hybrid_unif_gp,
                     x = x, 
                     method = "Nelder-Mead", 
                     control = list(maxit = 1000, fnscale = fscale)
                     )
  res = list()
  res$estimates = opt_result$par
  names(res$estimates) = c("U:lower bound","U:upper bound",
                           "GP:shape","GP:scale",
                           "mixing parameter")
  res$neg_LL = opt_result$value
  return(res)
}
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
q_hybrid_unif_gp = function(alpha, a, b, sh, sc, mix_par){
  # function that computes the alpha-quantile of a 
  # hybrid uniform-GP distribution with parameters
  # a, b, sh ,sc, mix_par
  
  if(alpha <= mix_par){
    return(qunif(alpha / mix_par, min = a, max = b))
  }else{
    return(mev::qgp((alpha - mix_par) / (1 - mix_par),
                    loc = b, scale = sc, shape = sh))
  }
}
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
q_hybrid_unif_weibull = function(alpha,a,b,sh,sc,mix_par){
  # function that computes the alpha-quantile of a 
  # hybrid uniform-GP distribution with parameters
  # a, b, sh ,sc, mix_par
  
  if(alpha <= mix_par){
    return(qunif(alpha / mix_par, min = a, max = b))
  }else{
    return(b + qweibull((alpha - mix_par) / (1 - mix_par), 
                        shape = sh, scale = sc))
    }
}
# ------------------------------------------------------------------------------

# ==============================================================================


# ==============================================================================
# Section "Multiple-threshold GP model"
# Thresholds defined as sample quantiles of a vector of observations x and 
# vector of probabilities q = (q_1, ..., q_m) defining thresholds u_1, ..., u_m;
# thresholds relative to u_1 defined as v_i = u_i - u_1, and distance between
# thresholds is w_i = u_{i+1} - u_i = v_{i+1} - v_i

# ------------------------------------------------------------------------------
dgpd_multiple_thresholds = function(y, x, q, theta){
  # function computing the density of the multiple-threshold GP distribution 
  # with parameter vector theta = (sigma1, xi_1, ..., xi_m)
  # y: excesses of threshold u_1
  
  quantiles = quantile(x, probs = q) 
  u = as.numeric(quantiles) 
  m = length(u) 
  v = u[1:m] - u[1]   
  w = c(diff(v), NA)
  
  y=y[y > u[1]]
  
  sigma1 = theta[1]                                        
  if (sigma1 <= 0) return(10^10)                                                # Need sigma_1 > 0
  xi = theta[-1]                                                                # (xi_1, ..., xi_m)
  sigma = sigma1 + cumsum(c(0, xi[-m] * w[-m]))                                 # (sigma_1, ..., sigma_m)
  phi = xi / sigma                                                              # (phi_1, ..., phi_m)
  if (any(1 + phi[-m] * w[-m] <= 0)) return(10^10)                              # Need all elements of 1 + phi * w / sigma > 0
  Ij = unlist(lapply(y, function(r) sum(r - v > 0)))                            # Interval indicators
  if (any(1 + phi[Ij] * (y - v[Ij]) <= 0)) return(10^10)                        # Need all elements of 1 + phi[Ij] * (y - v[Ij]) > 0
  aj = c(0, cumsum(log(1 + phi[-m] * w[-m]) / sigma[-m] / phi[-m]))             # - log(p_j), j = 1, ..., m
  pj = exp(-aj)                                                                 # P(Y > v_j), j = 1, ..., m
  bj = log(sigma)
  dj = log(1 + phi[Ij] * (y - v[Ij]))
  ej = log(1 + phi[Ij] * (y - v[Ij])) / sigma[Ij] / phi[Ij]
  exp(- (aj[Ij] + bj[Ij] + dj + ej))
}
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
qgpd_multiple_thresholds = function(x,q,alpha,theta){
  # function computing the alpha-quantile of the multiple-threshold 
  # GP distribution with parameter vector theta = (sigma1, xi_1, ..., xi_m)
  # for a probability alpha (scalar)
  
  quantiles = quantile(x, probs = q) 
  u = as.numeric(quantiles) 
  m = length(u) 
  y = x[x > u[1]] - u[1] 
  v = u[1:m] - u[1]    
  w = c(diff(v), NA)
  
  sigma1 = theta[1]                                                             # sigma_1
  if (sigma1 <= 0) return(10^10)                                                # Need sigma_1 > 0
  xi = theta[-1]                                                                # (xi_1, ..., xi_m)
  sigma = sigma1 + cumsum(c(0, xi[-m] * w[-m]))                                 # (sigma_1, ..., sigma_m)
  phi = xi / sigma                                                              # (phi_1, ..., phi_m)
  if (any(1 + phi[-m] * w[-m] <= 0)) return(10^10)                              # Need all elements of 1 + phi * w / sigma > 0
  Ij = unlist(lapply(y, function(r) sum(r - v > 0)))                            # Interval indicators
  if (any(1 + phi[Ij] * (y - v[Ij]) <= 0)) return(10^10)                        # Need all elements of 1 + phi[Ij] * (y - v[Ij]) > 0
  aj = c(0, cumsum(log(1 + phi[-m] * w[-m]) / sigma[-m] / phi[-m]))             # -log(p_j), j = 1, ..., m
  pj = unname(exp(-aj))                                                         # P(Y > v_j), j = 1, ..., m
  k = sum(1 - alpha - pj < 0)
  q = ifelse(xi[k] != 0,
             u[k] + (sigma[k] / xi[k]) * (((1 - alpha) / (pj[k]))^(-xi[k]) - 1),
             u[k] + sigma[k] * (log(1 - alpha) - log(pj[k])))
  q = unname(q)
  q
}
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
neg_log_likelihood = function(theta,x,q){
  # function computing the negative log-likelihood for the multiple-threshold 
  # GP model with parameter vector theta
  
  quantiles = quantile(x, probs = q) 
  u = as.numeric(quantiles) 
  m = length(u) 
  y = x[x > u[1]] - u[1] 
  v = u[1:m] - u[1]    
  w = c(diff(v), NA)
  
  sigma1 = theta[1]                                                             # sigma_1
  if (sigma1 <= 0) return(10^10)                                                # Need sigma_1 > 0
  xi = theta[-1]                                                                # (xi_1, ..., xi_m)
  sigma = sigma1 + cumsum(c(0, xi[-m] * w[-m]))                                 # (sigma_1, ..., sigma_m)
  phi = xi / sigma                                                              # (phi_1, ..., phi_m)
  if (any(1 + phi[-m] * w[-m] <= 0)) return(10^10)                              # Need all elements of 1 + phi * w / sigma > 0
  Ij = unlist(lapply(y, function(r) sum(r - v >0)))                             # Interval indicators
  if (any(1 + phi[Ij] * (y - v[Ij]) <= 0)) return(10^10)                        # Need all elements of 1 + phi[Ij] * (y - v[Ij]) > 0
  aj = c(0, cumsum(log(1 + phi[-m] * w[-m]) / sigma[-m] / phi[-m]))             # -log(p_j), j = 1, ..., m
  pj = exp(-aj)                                                                 # P(Y > v_j), j = 1, ..., m
  bj = log(sigma)
  dj = log(1 + phi[Ij] * (y - v[Ij]))
  ej = log(1 + phi[Ij] * (y - v[Ij])) / sigma[Ij] / phi[Ij]
  sum(aj[Ij] + bj[Ij] + dj + ej)
}
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
penalization = function(lambda, xi, type = "ridge"){
  # function computing the penalty on the first differences 
  # of the shape parameters
  
  m = length(xi)
  if(type == "ridge"){
      D = -diff(diag(m), differences = 1)
      0.5 * lambda * 
        t(D %*% matrix(xi, nrow = m, ncol = 1)) %*% 
        (D %*% matrix(xi, nrow = m, ncol = 1))
      }else if(type == "lasso"){
        lambda * sum(abs(D %*% matrix(xi, nrow = m, ncol = 1)))
        }else{stop("Only 'ridge' and 'lasso' penalization is available.")}
}
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
weighted_penalization = function(lambda, xi, theta_init, x, q){
  # function computing the weighted penalty on the first differences
  # of the shape parameters (ridge penalty only)
  
  quantiles = quantile(x, probs = q) 
  u = as.numeric(quantiles) 
  m = length(u) 
  y = x[x > u[1]] - u[1]
  v = u[1:m] - u[1]    
  w = c(diff(v), NA)
  
  D = -diff(diag(m))
  D2 = rbind(c(rep(0, ncol(D))), D)
  D2 = cbind(c(rep(0, nrow(D2))), D2)
  
  w_k = unname(1 / dgpd_multiple_thresholds(u, x, q, theta_init))
  W = diag(w_k)

  S_w = t(D2)%*%W%*%D2
  
  0.5 * lambda * 
    t(matrix(c(0, xi), nrow = m + 1, ncol = 1)) %*% 
    S_w %*% matrix(c(0, xi), nrow = m + 1, ncol = 1)
}
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
penalized_neg_log_likelihood = function(theta,lambda,x,q){
  xi = theta[-1]
  neg_log_likelihood(theta, x, q) + penalization(lambda, xi)
}
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
weighted_penalized_neg_log_likelihood = function(theta, lambda, theta_init, x, q){
  xi=theta[-1]
  neg_log_likelihood(theta, x, q) + 
    weighted_penalization(lambda, xi, theta_init, x, q)
}
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
exp.info.algebraic = function(x,z,q){
  quantiles = quantile(z, probs=q) 
  u = as.numeric(quantiles) 
  m = length(u) 
  y = z[z>u[1]]-u[1] # excesses of u[1]
  v = u[1:m]-u[1]    # thresholds relative to u[1] (to add m+1 index later..)
  w = c(diff(v),NA)
  
  # x : parameter vector: (sigma1,phi_1,...,phi_m)
  # v : thresholds relative to lowest threshold
  # w : differences between thresholds (w[m] not used)
  # m : number of thresholds
  
  sigma1 <- x[1]                                  # sigma_1
  phi <- x[-1]                                    # (phi_1,...,phi_m)
  sigma <- sigma1*c(1,cumprod(1+phi[-m]*w[-m]))   # (sigma_1,...,sigma_m)
  xi <- sigma*phi                                 # (xi_1,...,xi_m)
  aj <- c(0,cumsum(log(1+phi[-m]*w[-m])/sigma[-m]/phi[-m])) # -log(p_j), j=1,...,m
  pj <- exp(-aj)                                  # P(Y > v_j), j=1,...,m
  qj <- c(pj[-m]*(1-(1+phi[-m]*w[-m])^(-1/xi[-m])),pj[m]) # P(v_j < Y < v_{j+1}), j=1,...,m
  h <- phi*sigma/sigma1                           # (h_1,...,h_m)
  #
  ##################### Various integrals ...
  #
  I0b <- function(b,j){
    t1 <- (1+b*xi[j])^(-1)
    if (j==m) return(t1)
    t1-(1+phi[j]*w[j])^(-b-1/xi[j])*t1
  }
  #
  I1b <- function(b,j){
    t1 <- sigma[j]/(1+(b-1)*xi[j])/(1+b*xi[j])
    if (j==m) return(t1)
    t2 <- (1+phi[j]*w[j])^(-b-1/xi[j])/(1+b*xi[j])
    t3 <- (1+phi[j]*w[j])^(1-b-1/xi[j])/(1+(b-1)*xi[j])
    t1+(t2-t3)/phi[j]
  }
  #
  I2b <- function(b,j){
    t1 <- 2*sigma[j]^2/(1+(b-2)*xi[j])/(1+(b-1)*xi[j])/(1+b*xi[j])
    if (j==m) return(t1)
    t2 <- (1+phi[j]*w[j])^(2-b-1/xi[j])/(1+(b-2)*xi[j])
    t3 <- (1+phi[j]*w[j])^(1-b-1/xi[j])/(1+(b-1)*xi[j])
    t4 <- (1+phi[j]*w[j])^(-b-1/xi[j])/(1+b*xi[j])
    t1-(t2-2*t3+t4)/phi[j]^2
  }
  # 
  J <- function(j){
    t1 <- xi[j]
    if (j==m) return(t1)
    t2 <- -(1+phi[j]*w[j])^(-1/xi[j])*log(1+phi[j]*w[j])
    t3 <- -xi[j]*(1+phi[j]*w[j])^(-1/xi[j])   
    t1+t2+t3
  } 
  #
  # Derivatives of bj (constant w.r.t. y) ...
  #
  # Note: all contributions are zero apart from the diagonal elements
  #       for sigma1, phi_1,...,phi_{m-1}
  #
  b.exp.info <- matrix(0,m+1,m+1)              # matrix for expected information from b
  db.phi.phi <- c(unlist(lapply(1:(m-1),function(x)-w[x]^2*(1+phi[x]*w[x])^(-2)*sum(qj[(x+1):m]))),0)
  diag(b.exp.info) <- c(-1/sigma1^2,db.phi.phi)
  #
  # Derivatives of dj ...
  #
  # Note: all contributions are zero apart from the diagonal elements
  #       for phi_1,...,phi_m.
  #
  d.exp.info <- matrix(0,m+1,m+1)              # matrix for expected information from b
  dd.phi.phi <- unlist(lapply(1:m,function(x)-pj[x]*I2b(2,x)))
  diag(d.exp.info) <- c(0,dd.phi.phi)
  #
  # Derivatives of ej ...
  #
  e.exp.info <- matrix(0,m+1,m+1)              # matrix for expected information from b
  de.s1.s1 <- 2*sigma1^(-3)*sum(h^(-1)*pj*unlist(lapply(1:m,J)))
  #
  de.phi.fn <- function(x){
    temp1 <- h[x]^(-1)*( -phi[x]^(-1)*J(x)+I1b(1,x) )
    if (x==m) return(sigma1^(-1)*temp1*pj[x])
    temp2 <- -w[x]*(1+phi[x]*w[x])^(-1)*h[(x+1):m]^(-1)*unlist(lapply((x+1):m,J))
    sigma1^(-1)*(temp1*pj[x]+sum(temp2*pj[(x+1):m]))
  }
  #
  de.phi <- unlist(lapply(1:m,de.phi.fn))
  #
  de.phi.phi.fn <- function(x){
    temp1 <- h[x]^(-1)*( 2*phi[x]^(-2)*J(x)-2*phi[x]^(-1)*I1b(1,x)-I2b(2,x) )
    if (x==m) return(sigma1^(-1)*temp1*pj[x])
    temp2 <- 2*w[x]^2*(1+phi[x]*w[x])^(-2)*h[(x+1):m]^(-1)*unlist(lapply((x+1):m,J))  
    sigma1^(-1)*(temp1*pj[x]+sum(temp2*pj[(x+1):m]))
  }
  #
  de.phi.phi <- unlist(lapply(1:m,de.phi.phi.fn))
  #
  de.phil.phik.fn <- function(x){
    k <- x[1]; l <- x[2]
    temp1 <- h[k]^(-1)*w[l]*(1+w[l]*phi[l])^(-1)*(phi[k]^(-1)*J(k)-I1b(1,k))
    if (k==m) return(sigma1^(-1)*temp1*pj[k])
    temp2 <- w[k]*(1+phi[k]*w[k])^(-1)*w[l]*(1+phi[l]*w[l])^(-1)*h[(k+1):m]^(-1)*unlist(lapply((k+1):m,J))
    sigma1^(-1)*(temp1*pj[k]+sum(temp2*pj[(k+1):m]))
  }
  #
  k.vals <- rep(2:m,times=2:m-1)                                # k > l 
  l.vals <- unlist(lapply(2:m-1,function(x) seq(from=1,to=x)))
  kl.vals <- cbind(k.vals,l.vals)
  rev.kl.vals <- kl.vals[,2:1]
  if (m==2){
    kl.vals <- matrix(kl.vals,nrow=1,ncol=2) # if m=2 make kl.vals a matrix
    rev.kl.vals <- matrix(kl.vals[,2:1],nrow=1,ncol=2) # if m=2 make rev.kl.vals a matrix
  }
  de.phil.phik <- unlist(apply(kl.vals,1,de.phil.phik.fn))
  #
  diag(e.exp.info) <- c(de.s1.s1,de.phi.phi)
  e.exp.info[1,2:(m+1)] <- e.exp.info[2:(m+1),1] <- -sigma1^(-1)*de.phi
  e.exp.info[kl.vals+1] <- e.exp.info[rev.kl.vals+1] <- de.phil.phik
  #
  ##################### Derivatives of aj ...
  #
  a.exp.info <- matrix(0,m+1,m+1)              # matrix for expected information from b
  aj <- c(0,cumsum(log(1+phi[-m]*w[-m])/sigma[-m]/phi[-m])) # -log(p_j), j=1,...,m
  da.s1.s1 <- 2*sigma1^(-2)*sum(aj*qj)
  #
  da.phi.fn <- function(x){
    temp1 <- h[x]^(-1)*( -phi[x]^(-1)*log(1+phi[x]*w[x])+w[x]*(1+phi[x]*w[x])^(-1) ) # B(x,x)
    my.zeros <- numeric(m-1)
    my.zeros[x] <- temp1
    temp1 <- my.zeros
    if (x==(m-1)) return(sigma1^(-1)*sum(temp1*qj[-1]))
    if (x<(m-1)){
      ind <- (x+1):(m-1)
      temp2 <- -w[x]*(1+phi[x]*w[x])^(-1)*h[ind]^(-1)*log(1+phi[ind]*w[ind])
    }
    temp1[(x+1):(m-1)] <- temp2
    sigma1^(-1)*sum(cumsum(temp1)*qj[-1])
  }
  #
  da.phi <- c(unlist(lapply(1:(m-1),da.phi.fn)),0)
  da.s1.phi <- -sigma1^(-1)*da.phi
  #
  da.phi.phi.fn <- function(x){
    y.1v <- 1+phi[x]*w[x]
    y.v <- w[x]
    temp1 <- h[x]^(-1)*( 2*phi[x]^(-2)*log(y.1v)-2*phi[x]^(-1)*y.v*y.1v^(-1)-y.v^2*y.1v^(-2) )
    my.zeros <- numeric(m-1)
    my.zeros[x] <- temp1
    temp1 <- my.zeros
    if (x==(m-1)) return(sigma1^(-1)*sum(temp1*qj[-1]))
    if (x<(m-1)){
      ind <- (x+1):(m-1)
      temp2 <- 2*w[x]^2*(1+phi[x]*w[x])^(-2)*h[ind]^(-1)*log(1+phi[ind]*w[ind])
    }
    temp1[(x+1):(m-1)] <- temp2
    sigma1^(-1)*sum(cumsum(temp1)*qj[-1])
  }
  #
  da.phi.phi <- c(unlist(lapply(1:(m-1),da.phi.phi.fn)),0)
  #
  da.phil.phik.fn <- function(x){
    k <- x[1]; l <- x[2]
    y.1v <- 1+phi[k]*w[k]
    y.v <- w[k]
    temp1 <- h[k]^(-1)*w[l]*(1+w[l]*phi[l])^(-1)*(phi[k]^(-1)*log(y.1v)-y.v*y.1v^(-1))
    my.zeros <- numeric(m-1)
    my.zeros[k] <- temp1
    temp1 <- my.zeros
    if (k==(m-1)) return(sigma1^(-1)*sum(temp1*qj[-1]))
    if (k<(m-1)){
      ind <- (k+1):(m-1)
      temp2 <- w[k]*(1+phi[k]*w[k])^(-1)*w[l]*(1+phi[l]*w[l])^(-1)*h[ind]^(-1)*log(1+phi[ind]*w[ind])
    }
    temp1[(k+1):(m-1)] <- temp2
    sigma1^(-1)*sum(cumsum(temp1)*qj[-1])
  }
  #
  diag(a.exp.info) <- c(da.s1.s1,da.phi.phi)
  a.exp.info[1,2:(m+1)] <- a.exp.info[2:(m+1),1] <- da.s1.phi
  # 
  if (m>2){
    k.vals <- rep(2:(m-1),times=1:(m-2))                                # k > l 
    l.vals <- unlist(lapply(1:(m-2),function(x) seq(from=1,to=x)))
    kl.vals <- cbind(k.vals,l.vals)
    rev.kl.vals <- kl.vals[,2:1]
    if (m==3){
      kl.vals <- matrix(kl.vals,nrow=1,ncol=2) # if m=3 make kl.vals a matrix
      rev.kl.vals <- matrix(kl.vals[,2:1],nrow=1,ncol=2) # if m=3 make rev.kl.vals a matrix
    }
    da.phil.phik <- unlist(apply(kl.vals,1,da.phil.phik.fn))
    a.exp.info[kl.vals+1] <- a.exp.info[rev.kl.vals+1] <- da.phil.phik
  }
  #
  # Return observed information ...
  exp.info <- a.exp.info+b.exp.info+d.exp.info+e.exp.info
  #
  exp.info
}
# ------------------------------------------------------------------------------

# ==============================================================================


# ==============================================================================
# Section: Monotonic shape parameters: re-parametrized functions (argument gamma)

# ------------------------------------------------------------------------------
mono_neg_log_likelihood = function(gamma, S_mono, x, q){
  # function computing the negative log-likelihood of the multiple-threshold GP
  # model with monotonic shape parameters 
  
  # gamma: (sigma1, beta_1,..., beta_m) = 
  #        (sigma1, +-(xi_2 - xi_1), ..., +-(xi_m - xi_{m-1}))
  
  # gamma_exp: (sigma1, exp(beta_1),...,exp(beta_m))

  quantiles = quantile(x, probs = q) 
  u = as.numeric(quantiles) 
  m = length(u) 
  y = x[x > u[1]] - u[1] 
  v = u[1:m] - u[1]    
  w = c(diff(v), NA) 
  
  gamma_exp = c(gamma[1:2], exp(gamma[3:(m+1)]))
  theta = S_mono %*% gamma_exp
  sigma1 = theta[1]                                                             # sigma_1
  if (sigma1 <= 0) return(10^10)                                                # Need sigma_1 > 0
  xi = theta[-1]                                                                # (xi_1,...,xi_m)
  sigma = sigma1 + cumsum(c(0, xi[-m] * w[-m]))                                 # (sigma_1,...,sigma_m)
  phi = xi / sigma                                                              # (phi_1,...,phi_m)
  if (any(1 + phi[-m] * w[-m] <= 0)) return(10^10)                              # Need all elements of 1+phi*w/sigma > 0
  Ij = unlist(lapply(y, function(r) sum(r - v > 0)))                            # interval indicators
  if (any(1 + phi[Ij] * (y - v[Ij]) <= 0)) return(10^10)                        # Need all elements of 1+phi[Ij]*(y-v[Ij]) > 0
  aj = c(0, cumsum(log(1 + phi[-m] * w[-m]) / sigma[-m] / phi[-m]))             # -log(p_j), j=1,...,m
  pj = exp(-aj)                                                                 # P(Y > v_j), j=1,...,m
  bj = log(sigma)
  dj = log(1 + phi[Ij] * (y - v[Ij]))
  ej = log(1 + phi[Ij] * (y - v[Ij])) / sigma[Ij] / phi[Ij]
  sum(aj[Ij] + bj[Ij] + dj + ej)
}
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
mono_penalization = function(lambda, gamma, S_mono, q, type = "ridge"){
  m = length(q)
  Ip = diag(c(0,0,rep(1,m-1)))
  gamma_exp = c(gamma[1:2], exp(gamma[3:(m+1)]))
  if(type == "ridge"){
    0.5 * lambda * t(gamma_exp) %*% Ip %*% gamma_exp
    }else if(type == "lasso"){
      lambda * sum(Ip %*% gamma_exp)
      } else{stop("Only 'ridge' and 'lasso' penalities are available.")}
}
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
mono_penalized_neg_log_likelihood = function(gamma, S_mono, lambda, x, q){
  mono_neg_log_likelihood(gamma, S_mono, x, q) + 
    mono_penalization(lambda, gamma, S_mono, q)
}
# ------------------------------------------------------------------------------

# ==============================================================================


# ==============================================================================
# Section: simulation study

optimize_neg_log_likelihood = function(neg_LL_function, 
                                       x,
                                       q = c(seq(from = 0, to = 0.9, by = 0.3)), 
                                       lambda,
                                       theta_opt = rep(0, length(q) + 1), 
                                       gamma_opt = rep(0, length(q) + 1),
                                       theta_init = rep(0, length(q) + 1), 
                                       method = "unpenalized"){
  
  if(method == "unpenalized"){
    fscale = neg_LL_function(theta_init, x, q)
    tryCatch({
      result = optim(par = theta_init, neg_LL_function, x = x, q = q,
                     hessian = F, method = "BFGS", 
                     control = list(maxit = 100, fnscale = fscale))
      }, error = function(e){
        result = optim(par = theta_init, neg_LL_function, x = x, q = q,
                       hessian = F, method = "CG", 
                       control = list(maxit = 100, fnscale = fscale)) 
    })
  }
  
  if(method == "penalized"){
    fscale = neg_LL_function(theta = theta_opt, lambda = lambda, x = x, q = q)
    tryCatch({
      result = optim(par = theta_opt, fn = neg_LL_function,
                     lambda = lambda, x = x, q = q,
                     hessian = F, method = "BFGS", 
                     control = list(maxit = 100, fnscale = fscale))
      }, error = function(e){
        result = optim(par = theta_opt, fn = neg_LL_function,
                       lambda = lambda, x = x, q = q,
                       hessian = F, method = "CG", 
                       control = list(maxit = 100, fnscale = fscale))
    })
  }
  
  if(method == "penalized_weighted"){
    fscale =  neg_LL_function(theta = theta_opt, lambda = lambda, 
                              theta_init = theta_init, x = x, q = q)
    tryCatch({
      result = optim(par = theta_opt, fn = neg_LL_function,
                     lambda = lambda, theta_init = theta_init, x = x, q = q,
                     hessian = F, method = "BFGS", 
                     control = list(maxit = 100, fnscale = fscale))
    }, error = function(e){
      result = optim(par = theta_opt, fn = neg_LL_function,
                     lambda = lambda, theta_init = theta_init, x = x, q = q,
                     hessian = F, method = "CG", 
                     control = list(maxit = 100, fnscale = fscale))
    })
  }
  
  if(method == "penalized_increasing"){
    m = length(q)
    matrix_of_ones = matrix(1, nrow = m, ncol = m)
    S_inc = matrix_of_ones * lower.tri(matrix_of_ones, diag = TRUE)
    S_inc = rbind(c(rep(0, m)), S_inc)
    S_inc = cbind(c(1, rep(0, m)), S_inc)
    
    fscale = neg_LL_function(gamma = gamma_opt, S_mono = S_inc, 
                             lambda = lambda, x = x, q = q)
    tryCatch({
      result = optim(par = gamma_opt, fn = neg_LL_function,
                     S_mono = S_inc,
                     lambda = lambda, x = x, q = q,
                     hessian = F, method = "BFGS", 
                     control = list(maxit = 100, fnscale = fscale))
    }, error=function(e){
      result = optim(par = gamma_opt, fn = neg_LL_function,
                     S_mono = S_inc,
                     lambda = lambda, x = x, q = q,
                     hessian = F, method = "CG", 
                     control = list(maxit = 100, fnscale = fscale))
    })
    gamma_hat_lambda = result$par
    gamma_exp_hat_lambda = c(gamma_hat_lambda[1:2], exp(gamma_hat_lambda[3:(m+1)]))
    theta_hat_lambda = S_inc %*% gamma_exp_hat_lambda
    result$par = theta_hat_lambda
  }
  
  if(method == "penalized_decreasing"){
    m = length(q)
    matrix_of_ones = matrix(1, nrow = m, ncol = m)
    S_dec = - matrix_of_ones * lower.tri(matrix_of_ones, diag = TRUE)
    S_dec[, 1] = rep(1, m)
    S_dec = rbind(c(rep(0, m)), S_dec)
    S_dec = cbind(c(1, rep(0, m)), S_dec)
    
    fscale = neg_LL_function(gamma = gamma_opt, S_mono = S_dec,
                             lambda = lambda, x = x, q = q)
    tryCatch({
      result = optim(par = gamma_opt, fn = neg_LL_function,
                     S_mono = S_dec,
                     lambda = lambda, x = x, q = q,
                     hessian = F, method = "BFGS", 
                     control = list(maxit = 100, fnscale = fscale))
    }, error = function(e){
      result = optim(par = gamma_opt, fn = neg_LL_function,
                     S_mono = S_dec,
                     lambda = lambda, x = x, q = q,
                     hessian = F, method = "CG", 
                     control = list(maxit = 100, fnscale = fscale))
    })
    
    gamma_hat_lambda = result$par
    gamma_exp_hat_lambda = c(gamma_hat_lambda[1:2], exp(gamma_hat_lambda[3:(m+1)]))
    theta_hat_lambda = S_dec %*% gamma_exp_hat_lambda
    result$par = theta_hat_lambda
  }
  
  if(method != "unpenalized" && method != "penalized" && 
     method != "penalized_weighted" && method != "penalized_increasing" &&
     method != "penalized_decreasing"){
    stop("Invalid method")
  } else{
      return(result)
    }
  
}

compute_marg_neg_LL = function(theta_hat_lambda, x, q, lambda, S, optim_result){
  # function that computes the value of the log-marginal negative likelihood
  
  quantiles = quantile(x, probs = q) 
  u = as.numeric(quantiles) 
  m = length(u) 
  y = x[x > u[1]] - u[1] 
  v = u[1:m] - u[1]    
  w = c(diff(v), NA)  
  
  xi_hat = theta_hat_lambda[-1]                                                                
  sigma_hat = theta_hat_lambda[1] + cumsum(c(0, xi_hat[-m] * w[-m]))            
  phi_hat = xi_hat / sigma_hat
  psi_hat = c(sigma_hat[1], phi_hat)
  
  exp_info = length(y) * exp.info.algebraic(psi_hat, x, q)
  
  svd = svd(exp_info + lambda * S)
  eigen_H = eigen(exp_info + lambda * S)
  
  if(prod(svd$d)>.Machine$double.xmin){
    return(optim_result$value + 0.5 * log(prod(svd$d)) - 0.5 * (m - 1) * log(lambda))
  }else if(det(exp_info + lambda * S)>.Machine$double.xmin){
    return(optim_result$value + 0.5 * log(det(exp_info + lambda * S)) - 0.5 * (m-1) * log(lambda))
  }else if(prod(eigen_H$values)>.Machine$double.xmin){
    return(optim_result$value + 0.5 * log(prod(eigen_H$values)) - 0.5 * (m-1) * log(lambda))
  }else{
    return(.Machine$double.xmax)
  }
}

# ------------------------------------------------------------------------------
marg_neg_LL = function(log_lambda, x, q, 
                       theta_opt = rep(0, length(q) + 1), 
                       gamma_opt = rep(0, length(q) + 1), 
                       method = "penalized"){
  
  # method = "penalized" / "penalized_weighted" / "penalized_increasing" / "penalized_decreasing"
  
  lambda=10^(log_lambda)
  
  # ---- Data-related quantities ----
  quantiles = quantile(x, probs = q) 
  u = as.numeric(quantiles) 
  m = length(u) 
  y = x[x > u[1]] - u[1] 
  v = u[1:m] - u[1]    
  w = c(diff(v), NA)  
  n = length(y)
  
  # ---- Initialize penalty matrix ----
  D = -diff(diag(m), differences = 1)
  D2 = rbind(c(rep(0, ncol(D))), D)
  D2 = cbind(c(rep(0, nrow(D2))), D2)
  S = crossprod(D2)                                                             
  
  # ---- Monotonicity: transformation matrices (gamma vs theta) ----
  matrix_of_ones = matrix(1, nrow = m, ncol = m)

  # S_inc and S_dec are (m+1)x(m+1)
  S_inc = matrix_of_ones * lower.tri(matrix_of_ones, diag = TRUE)
  S_inc = rbind(c(rep(0, m)), S_inc)
  S_inc = cbind(c(1, rep(0, m)), S_inc)

  S_dec = - matrix_of_ones * lower.tri(matrix_of_ones, diag = TRUE)
  S_dec[, 1] = rep(1, m)
  S_dec = rbind(c(rep(0, m)), S_dec)
  S_dec = cbind(c(1, rep(0, m)), S_dec)

  # ---- Case by case ----
  if(method == "penalized"){
    # ---- Optimization of log-likelihood ----
    result = optimize_neg_log_likelihood(penalized_neg_log_likelihood, 
                                         theta_opt = theta_opt, 
                                         lambda = lambda, x = x, q = q,
                                         method = "penalized")
    
    theta_hat_lambda = result$par
    
    return(compute_marg_neg_LL(theta_hat_lambda = theta_hat_lambda,
                               x = x, q = q, lambda = lambda, S = S,
                               optim_result = result))
    
  } 
  else if(method == "penalized_weighted"){
    
    # ---- Initial fit ----
    initial_fit = mev::fit.gpd(y)
    scale_init = initial_fit$estimate[1]
    shape_init = rep(initial_fit$estimate[2], m)
    theta_init = c(scale_init, shape_init)
    
    # ---- Weighted penalty matrix ----
    w_k = unname(1 / dgpd_multiple_thresholds(u, x, q, theta_init))
    W = diag(w_k)
    S_w = t(D2) %*% W %*% D2  
    
    # ---- Optimize neg. log-likelihood ----
    result = optimize_neg_log_likelihood(weighted_penalized_neg_log_likelihood,
                                         theta_opt = theta_opt, 
                                         theta_init = theta_init,
                                         x = x, q = q, lambda = lambda,
                                         method = "penalized_weighted")
    
    theta_hat_lambda = result$par
    
    return(compute_marg_neg_LL(theta_hat_lambda = theta_hat_lambda,
                               x = x, q = q, lambda = lambda, S = S_w,
                               optim_result = result))
  }
  else if(method == "penalized_increasing"){
    result = optimize_neg_log_likelihood(mono_penalized_neg_log_likelihood, 
                                         gamma_opt = gamma_opt,
                                         lambda = lambda, x = x, q = q,
                                         method = "penalized_increasing")
    gamma_hat_lambda = result$par
    gamma_exp_hat_lambda = c(gamma_hat_lambda[1:2], exp(gamma_hat_lambda[3:(m+1)]))
    theta_hat_lambda = S_inc %*% gamma_exp_hat_lambda
    result$par = theta_hat_lambda
    
    return(compute_marg_neg_LL(theta_hat_lambda = theta_hat_lambda,
                               x = x, q = q, lambda = lambda, S = S,
                               optim_result = result))
    
  } 
  else if(method == "penalized_decreasing"){
    
    result = optimize_neg_log_likelihood(mono_penalized_neg_log_likelihood, 
                                         gamma_opt = gamma_opt,
                                         lambda = lambda, x = x, q = q,
                                         method = "penalized_decreasing")
    
    gamma_hat_lambda = result$par
    gamma_exp_hat_lambda = c(gamma_hat_lambda[1:2], exp(gamma_hat_lambda[3:(m+1)]))
    theta_hat_lambda = S_dec %*% gamma_exp_hat_lambda
    result$par = theta_hat_lambda
    
    return(compute_marg_neg_LL(theta_hat_lambda = theta_hat_lambda,
                               x = x, q = q, lambda = lambda, S,
                               optim_result = result))
  } 
  else {stop("Invalid method value")}
}
# ------------------------------------------------------------------------------



# ------------------------------------------------------------------------------
m.t.GP.fit = function(x, q, 
                      unpenalized = TRUE, penalized = TRUE,
                      penalized_weighted = FALS,
                      penalized_increasing = FALSE, 
                      penalized_decreasing = FALSE, 
                      val = c(10^seq(-2, 16, by = 1))){
  
  # function that fits a multiple-threshold GP model (unpenalized / penalized /
  # penalized weighted / penalized increasing / penalized decreasing) to a 
  # unique dataset x 
  
  res = list() # list to stock results
  
  # ---- Data-related quantities ----
  quantiles = quantile(x, probs=q) 
  u = as.numeric(quantiles) 
  m = length(u) 
  y = x[x>u[1]]-u[1] 
  v = u[1:m]-u[1]    
  w = c(diff(v),NA)  
  n = length(y)
  
  # ---- Initialize penalty matrices ----
  
  D = -diff(diag(m), differences = 1)
  D2 = rbind(c(rep(0,ncol(D))),D)
  D2 = cbind(c(rep(0,nrow(D2))),D2)
  S = crossprod(D2)     # dim: (m+1)x(m+1)
  
  l = length(val)
  
  res$info$lambda_grid = val
  res$info$sample_size = length(x)
  res$info$nexc = n
  res$info$quantiles = q
  res$info$data = x
  
  # ---- Initial fit ----
  initial_fit = mev::fit.gpd(y)    
  scale_init = initial_fit$estimate[1]
  shape_init = rep(initial_fit$estimate[2], m)
  theta_init = c(scale_init, shape_init) 
  res$singleGP$estimates = initial_fit$estimate
  res$singleGP$std_err = initial_fit$std.err
  res$singleGP$neg_LL = initial_fit$nllh
  
  res$singleGP$AIC = 2 * initial_fit$nllh + 2
  res$singleGP$AICc = 2 * initial_fit$nllh + 2 + (2 * 2 * 3)/(n - 3)
  res$singleGP$BIC = 2 * initial_fit$nllh + 2 * log(n)
  
  if(unpenalized == TRUE){
    # ---- MLE without penalty ----
    result = optimize_neg_log_likelihood(neg_log_likelihood, x = x, q = q,
                                         theta_init = theta_init,
                                         method = "unpenalized")
    
    res$unpenalized$theta = result$par
    res$unpenalized$neg_LL = result$value
    
    xi_hat = result$par[-1]                                                               
    sigma_hat = result$par[1] + cumsum(c(0, xi_hat[-m] * w[-m]))                                 
    phi_hat = xi_hat / sigma_hat
    psi_hat = c(sigma_hat[1], phi_hat)
    
    exp_info = length(y) * exp.info.algebraic(psi_hat, x, q)
    svd = svd(exp_info)
    
    tryCatch({
      std_err_best = diag(solve(exp_info))
    }, error = function(e){
      U = svd$u
      L = diag(svd$d)
      V = svd$v
      L_inv = diag(1 / svd$d)
      inv_svd = V %*% L_inv %*% t(U)
      std_err_best = diag(inv_svd)
    })
    
    if(all(std_err_best >= 0)){
      std_err_best = sqrt(std_err_best)
    }else{
        std_err_best[which(std_err_best < 0)] = NA
        std_err_best[which(std_err_best >= 0)] = sqrt(std_err_best[which(std_err_best >= 0)])
    }
    
    res$unpenalized$exp_info = exp_info
    res$unpenalized$svd_exp_info = svd
    tryCatch({
      res$unpenalized$std_err = std_err_best
    }, error=function(e){
      res$unpenalized$std_err = NA
    })
    
    res$unpenalized$AIC = 2*temp$value+2*(m+1)
    res$unpenalized$AICc = 2*temp$value+2*(m+1)+(2*(m+1)*(m+2))/(n-m-2)
    res$unpenalized$BIC = 2*temp$value+(m+1)*log(n)
  }
  
  # ---- CASE 0: UNCONSTRAINED SHAPE PARAMETERS ----
  
  if(penalized == TRUE){
    # ---- Initialize vectors ----
    index = 0
    penalized_neg_LL = rep(NA,l)
    neg_LL = rep(NA,l)
    
    df = rep(NA,l)
    AIC = rep(NA,l)
    AICc = rep(NA,l)
    BIC = rep(NA,l) 
    marginal_neg_LL= rep(NA,l) 
    
    theta_hat_list = list() 
    exp_info_list = list()
    svd_list = list()
    std_err_list = list()
    theta_loop = theta_init
    
    # ---- "for" loop over lambda(s) ----
    for(lambda in val){
      
      index = index + 1
      
      # ---- Case 1: Standard penalty ----
      
      result = optimize_neg_log_likelihood(penalized_neg_log_likelihood,
                                           theta_opt = theta_loop,
                                           x = x, q = q, lambda = lambda,
                                           method = "penalized")
      
      theta_hat_lambda = result$par
      theta_hat_list[[index]] = theta_hat_lambda
      
      neg_LL[index] = neg_log_likelihood(theta_hat_lambda, x, q)
      penalized_neg_LL[index] = result$value
      
      xi_hat = theta_hat_lambda[-1]                                                                # (xi_1,...,xi_m)
      sigma_hat = theta_hat_lambda[1] + cumsum(c(0, xi_hat[-m] * w[-m]))                                      # (sigma_1,...,sigma_m)
      phi_hat = xi_hat/sigma_hat
      psi_hat = c(sigma_hat[1], phi_hat)
      
      exp_info = length(y) * exp.info.algebraic(psi_hat, x, q)
      svd = svd(exp_info + lambda * S)
      
      tryCatch({
        prod = solve(exp_info + lambda * S) %*% exp_info
        std_err = diag(solve(exp_info + lambda * S))
      }, error = function(e){
        U = svd$u
        L = diag(svd$d)
        V = svd$v
        L_inv = diag(1/svd$d)
        inv_svd = V %*% L_inv %*% t(U)
        prod = inv_svd %*% exp_info
        std_err = diag(inv_svd)
      })
      
      if(all(std_err >= 0)){
        std_err = sqrt(std_err)
      }else{
        std_err[which(std_err < 0)] = NA
        std_err[which(std_err >= 0)] = sqrt(std_err[which(std_err >= 0)])
      }
      
      exp_info_list[[index]] = exp_info
      svd_list[[index]] = svd
      std_err_list[[index]] = std_err
      
      eigen_H = eigen(exp_info + lambda * S)
      
      df[index] = sum(diag(prod))
      AIC[index] = 2 * neg_LL[index] + 2 * df[index]
      AICc[index] = AIC[index] + (2 * df[index] * (df[index] + 1))/(n - df[index] - 1)
      BIC[index] = 2 * neg_LL[index] + df[index] * log(n)
      
      if(prod(svd$d)>.Machine$double.xmin){
        marginal_neg_LL[index] = result$value + 0.5 * log(prod(svd$d)) - 0.5 * (m - 1) * log(lambda)
      }else if(det(exp_info+lambda*S)>.Machine$double.xmin){
        marginal_neg_LL[index] = result$value + 0.5 * log(det(exp_info + lambda * S)) - 0.5 * (m - 1) * log(lambda)
      }else if(prod(eigen_H$values)>.Machine$double.xmin){
        marginal_neg_LL[index] = result$value + 0.5 * log(prod(eigen_H$values)) - 0.5 * (m - 1) * log(lambda)
      }else{
        marginal_neg_LL[index] = .Machine$double.xmax
      }
      
      theta_loop = theta_hat_lambda
    }
    
    # ---- Optimal penalties ----
    I_opt = min(which(-marginal_neg_LL >= max(-marginal_neg_LL) - 0.01))
    #I_opt = which.max(-marginal_neg_LL)  
    theta_opt = theta_hat_list[[I_opt]]
    
    lower = log10(val[max(I_opt - 2, 1)])
    upper = log10(val[min(I_opt + 2, l)])
    marg_neg_LL_optimized = optimize(marg_neg_LL,
                                     interval = c(lower, upper), 
                                     tol = .Machine$double.eps^0.25,
                                     maximum = FALSE,
                                     method = "penalized",
                                     x = x, q = q, theta_opt = theta_opt)
    
    lambda_best = 10^(marg_neg_LL_optimized$minimum)
    
    result = optimize_neg_log_likelihood(penalized_neg_log_likelihood,
                                         theta_opt = theta_opt,
                                         x = x, q = q, lambda = lambda_best,
                                         method = "penalized")
    theta_best = result$par
    neg_LL_best = neg_log_likelihood(theta_best, x, q)
    penalized_neg_LL_best = result$value
    
    xi_hat = theta_best[-1]                                                      
    sigma_hat = theta_best[1] + cumsum(c(0, xi_hat[-m] * w[-m]))                  
    phi_hat = xi_hat/sigma_hat
    psi_hat = c(sigma_hat[1], phi_hat)
    
    exp_info = length(y) * exp.info.algebraic(psi_hat, x, q)
    svd = svd(exp_info + lambda_best * S)
    
    tryCatch({
      prod = solve(exp_info + lambda_best * S) %*% exp_info
      std_err_best = diag(solve(exp_info + lambda_best * S))
    }, error = function(e){
      
      U = svd$u
      L = diag(svd$d)
      V = svd$v
      L_inv = diag(1/svd$d)
      inv_svd = V %*% L_inv %*% t(U)
      prod = inv_svd %*% exp_info
      std_err_best = diag(inv_svd)
    })
    
    if(all(std_err_best >= 0)){
      std_err_best = sqrt(std_err_best)
    }else{
      std_err_best[which(std_err_best < 0)] = NA
      std_err_best[which(std_err_best >= 0)] = sqrt(std_err_best[which(std_err_best >= 0)])
    }
    
    res$penalized$exp_info = exp_info
    res$penalized$svd_exp_info = svd
    
    tryCatch({
      res$penalized$std_err = std_err_best
    }, error=function(e){
      res$penalized$std_err = NA
    })
    
    eigen_H = eigen(exp_info + lambda_best * S)
    
    df_best = sum(diag(prod))
    AIC_best = 2 * neg_LL_best + 2 * df_best
    AICc_best = AIC_best + (2 * df_best * (df_best + 1))/(n - df_best - 1)
    BIC_best = 2 * neg_LL_best + df_best * log(n)
    
    if(prod(svd$d)>.Machine$double.xmin){
      marginal_neg_LL_best = result$value + 0.5 * log(prod(svd$d)) - 0.5 * (m - 1) * log(lambda)
    }else if(det(exp_info+lambda_best*S)>.Machine$double.xmin){
      marginal_neg_LL_best = result$value + 0.5 * log(det(exp_info + lambda_best * S)) - 0.5 * (m - 1) * log(lambda_best)
    } else if(prod(eigen_H$values)>.Machine$double.xmin){
      marginal_neg_LL_best = result$value + 0.5 * log(prod(eigen_H$values)) - 0.5 * (m - 1) * log(lambda_best)
    }else{
      marginal_neg_LL_best = .Machine$double.xmax
    }
    
    
    # ---- Save optimal results ----
    res$penalized$theta = theta_best
    res$penalized$neg_LL = neg_LL_best
    res$penalized$penalized_neg_LL = penalized_neg_LL_best
    res$penalized$marginal_neg_LL = marginal_neg_LL_best
    res$penalized$AIC = AIC_best
    res$penalized$AICc = AICc_best
    res$penalized$BIC = BIC_best
    res$penalized$lambda = lambda_best
    res$penalized$df = df_best
    
    # ---- Save other results ----
    res$penalized$AIC_vec = AIC
    res$penalized$AICc_vec = AICc 
    res$penalized$BIC_vec = BIC
    res$penalized$df_vec = df
    res$penalized$marginal_neg_LL_vec = marginal_neg_LL
    res$penalized$neg_LL_vec = neg_LL 
    res$penalized$penalized_neg_LL_vec = penalized_neg_LL
    res$penalized$theta_hat_list = theta_hat_list 
    res$penalized$exp_info_list = exp_info_list
    res$penalized$svd_exp_info_list = svd_list
    res$penalized$std_err_list = std_err_list
  }
  
  # ---- CASE 0.1 : UNCONSTRAINED SHAPE PARAMETERS AND WEIGHTED PENALTY ----
  if(penalized_weighted == TRUE){
    
    # ---- Weighted penalty matrix ----
    w_k = unname(1/dgpd_multiple_thresholds(u,x,q,theta_init))
    W = diag(w_k)
    S_w = t(D2)%*%W%*%D2  # dim: (m+1)x(m+1)
    
    # ---- Initialize vectors ----
    index = 0
    w_penalized_neg_LL = rep(NA,l)
    w_neg_LL = rep(NA,l)
    
    w_df = rep(NA,l)
    w_AIC = rep(NA,l)
    w_AICc = rep(NA,l)
    w_BIC = rep(NA,l)
    w_marginal_neg_LL= rep(NA,l)
    
    w_theta_hat_list = list()
    w_exp_info_list = list()
    w_svd_list = list()
    w_std_err_list = list()
    w_theta_loop = theta_init 
    
    for(lambda in val){
      index = index + 1
      
      # ---- Case 2: Weighted penalty ---- 
      result = optimize_neg_log_likelihood(weighted_penalized_neg_log_likelihood,
                                           theta_opt = w_theta_loop,
                                           theta_init = theta_init, 
                                           x = x, q = q, lambda = lambda,
                                           method = "penalized_weighted")
      
      theta_hat_lambda = result$par
      w_theta_hat_list[[index]] = theta_hat_lambda 
      
      w_neg_LL[index] = neg_log_likelihood(theta_hat_lambda,x,q)
      w_penalized_neg_LL[index] = result$value
      
      xi_hat = theta_hat_lambda[-1]                                                                # (xi_1,...,xi_m)
      sigma_hat = theta_hat_lambda[1]+cumsum(c(0, xi_hat[-m] * w[-m]))                                      # (sigma_1,...,sigma_m)
      phi_hat = xi_hat/sigma_hat
      psi_hat = c(sigma_hat[1], phi_hat)
      
      exp_info = length(y) * exp.info.algebraic(psi_hat, x, q)
      svd = svd(exp_info + lambda * S_w)
      
      tryCatch({
        prod = solve(exp_info + lambda * S_w) %*% exp_info
        std_err = diag(solve(exp_info + lambda * S_w))
      }, error = function(e){
        U = svd$u
        L = diag(svd$d)
        V = svd$v
        L_inv = diag(1/svd$d)
        inv_svd = V %*% L_inv %*% t(U)
        prod = inv_svd %*% exp_info
        std_err = diag(inv_svd)
      })
      
      if(all(std_err >= 0)){
        std_err = sqrt(std_err)
      }else{
        std_err[which(std_err < 0)] = NA
        std_err[which(std_err >= 0)] = sqrt(std_err[which(std_err >= 0)])
      }
      
      w_exp_info_list[[index]] = exp_info
      w_svd_list[[index]] = svd
      w_std_err_list[[index]] = std_err
      
      eigen_H = eigen(exp_info + lambda * S_w)
      
      w_df[index] = sum(diag(prod))
      w_AIC[index] = 2 * w_neg_LL[index]+ 2 * w_df[index]
      w_AICc[index] = w_AIC[index] + (2 * w_df[index] * (w_df[index] + 1))/(n - w_df[index] - 1)
      w_BIC[index] = 2 * w_neg_LL[index] + w_df[index] * log(n) 
      
      if(prod(svd$d)>.Machine$double.xmin){
        w_marginal_neg_LL[index] = result$value + 0.5 * log(prod(svd$d)) - 0.5 * (m - 1) * log(lambda)
      }else if(det(exp_info+lambda*S_w)>.Machine$double.xmin){
        w_marginal_neg_LL[index] = result$value + 0.5 * log(det(exp_info + lambda * S_w)) - 0.5 * (m - 1) * log(lambda)
      }else if(prod(eigen_H$values)>.Machine$double.xmin){
        w_marginal_neg_LL[index] = result$value + 0.5 * log(prod(eigen_H$values)) - 0.5 * (m - 1) * log(lambda)
      }else{
        w_marginal_neg_LL[index] = .Machine$double.xmax
      }
      w_theta_loop = theta_hat_lambda
    }
    
    # ---- Optimal penalties ----
    w_I_opt = min(which(-w_marginal_neg_LL >= max(-w_marginal_neg_LL) - 0.01))
    #w_I_opt = which.max(-w_marginal_neg_LL)
    w_theta_opt = w_theta_hat_list[[w_I_opt]]
    
    lower = log10(val[max(w_I_opt - 2, 1)])
    upper = log10(val[min(w_I_opt + 2, l)])
    w_marg_neg_LL_optimized = optimize(marg_neg_LL,
                                       interval = c(lower, upper), 
                                       tol = .Machine$double.eps^0.25,
                                       maximum = FALSE,
                                       x = x, q = q, theta_opt = w_theta_opt, 
                                       method = "penalized_weighted")
    
    w_lambda_best = 10^(marg_neg_LL_optimized$minimum)
    
    result = optimize_neg_log_likelihood(weighted_penalized_neg_log_likelihood,
                                         theta_opt = w_theta_opt,
                                         theta_init = theta_init, 
                                         x = x, q = q, lambda = w_lambda_best,
                                         method = "penalized_weighted")
    
    w_theta_best = result$par
    w_neg_LL_best = neg_log_likelihood(w_theta_best,x,q)
    w_penalized_neg_LL_best = result$value
    
    xi_hat = w_theta_best[-1]                                                                # (xi_1,...,xi_m)
    sigma_hat = w_theta_best[1] + cumsum(c(0, xi_hat[-m] * w[-m]))                                      # (sigma_1,...,sigma_m)
    phi_hat = xi_hat/sigma_hat
    psi_hat = c(sigma_hat[1], phi_hat)
    
    exp_info = length(y) * exp.info.algebraic(psi_hat, x, q)
    svd = svd(exp_info + w_lambda_best * S_w)
    
    tryCatch({
      prod = solve(exp_info + w_lambda_best * S_w) %*% exp_info
      std_err_best = diag(solve(exp_info + w_lambda_best * S_w))
    }, error = function(e){
      U = svd$u
      L = diag(svd$d)
      V = svd$v
      L_inv = diag(1/svd$d)
      inv_svd = V %*% L_inv %*% t(U)
      prod = inv_svd %*% exp_info
      std_err_best = diag(inv_svd)
    })
    
    if(all(std_err_best >= 0)){
      std_err_best = sqrt(std_err_best)
    }else{
      std_err_best[which(std_err_best < 0)] = NA
      std_err_best[which(std_err_best >= 0)] = sqrt(std_err_best[which(std_err_best >= 0)])
    }
    
    res$weighted$exp_info = exp_info
    res$weighted$svd_exp_info = svd
    
    tryCatch({
      res$weighted$std_err = std_err_best
    }, error=function(e){
      res$weighted$std_err = NA
    })
    
    eigen_H = eigen(exp_info + w_lambda_best * S_w)
    
    w_df_best = sum(diag(prod))
    w_AIC_best = 2 * w_neg_LL_best + 2 * w_df_best
    w_AICc_best = w_AIC_best + (2 * w_df_best * (w_df_best + 1))/(n - w_df_best - 1)
    w_BIC_best = 2 * w_neg_LL_best + w_df_best * log(n)
    
    if(prod(svd$d)>.Machine$double.xmin){
      w_marginal_neg_LL_best = result$value + 0.5 * log(prod(svd$d)) - 0.5 * (m - 1) * log(lambda)
    }else if(det(exp_info+w_lambda_best*S_w)>.Machine$double.xmin){
      w_marginal_neg_LL_best = result$value + 0.5 * log(det(exp_info + w_lambda_best * S_w)) - 0.5 * (m - 1) * log(w_lambda_best)
    }else if(prod(eigen_H$values)>.Machine$double.xmin){
      w_marginal_neg_LL_best = result$value + 0.5 * log(prod(eigen_H$values)) - 0.5 * (m - 1) * log(w_lambda_best)
    }
    
    
    # ---- Save optimal results ----
    res$weighted$theta = w_theta_best
    res$weighted$neg_LL = w_neg_LL_best
    res$weighted$penalized_neg_LL = w_penalized_neg_LL_best
    res$weighted$marginal_neg_LL = w_marginal_neg_LL_best
    res$weighted$AIC = w_AIC_best
    res$weighted$AICc = w_AICc_best
    res$weighted$BIC = w_BIC_best
    res$weighted$lambda = w_lambda_best
    res$weighted$df = w_df_best
    
    # ---- Save other results ----
    res$weighted$AIC_vec = w_AIC
    res$weighted$AICc_vec = w_AICc 
    res$weighted$BIC_vec = w_BIC
    res$weighted$df_vec = w_df
    res$weighted$marginal_neg_LL_vec = w_marginal_neg_LL
    res$weighted$neg_LL_vec = w_neg_LL 
    res$weighted$penalized_neg_LL_vec = w_penalized_neg_LL
    res$weighted$theta_hat_list = w_theta_hat_list 
    res$weighted$exp_info_list = w_exp_info_list
    res$weighted$svd_exp_info_list = w_svd_list
    res$weighted$std_err_list = w_std_err_list
  }
  
  # ---- CASE 1: MONOTONICALLY INCREASING SHAPE PARAMETERS ----
  if(penalized_increasing == TRUE){
    
    matrix_of_ones = matrix(1, nrow = m, ncol = m)
    # S_inc and S_dec are (m+1)x(m+1)
    S_inc = matrix_of_ones * lower.tri(matrix_of_ones, diag = TRUE)
    S_inc = rbind(c(rep(0,m)),S_inc)
    S_inc = cbind(c(1,rep(0,m)),S_inc)
  
    Ip = diag(c(0,0,rep(1,m-1)))
    
    # ---- Initialize vectors ----
    index = 0
    inc_penalized_neg_LL = rep(NA,l) 
    inc_neg_LL = rep(NA,l) 
    
    inc_df = rep(NA,l) 
    inc_AIC = rep(NA,l)
    inc_AICc = rep(NA,l)
    inc_BIC = rep(NA,l) 
    inc_marginal_neg_LL= rep(NA,l)
    
    inc_theta_hat_list = list()
    inc_gamma_hat_list = list() 
    inc_exp_info_list = list()
    inc_svd_list = list()
    inc_std_err_list = list()
    
    sgn = 1 # sign + for increasing (- for decreasing)
    S_mono = S_inc
    gamma_init = c(theta_init[1:2], seq(from = log(0.2), to = log(3), length.out = m-1))
    gamma_loop = gamma_init
    
    # ---- "for" loop over lambda(s) ----
    for(lambda in val){
      index = index + 1
      
      result = optimize_neg_log_likelihood(mono_penalized_neg_log_likelihood,
                                           gamma_opt = gamma_loop,
                                           x = x, q = q, lambda = lambda,
                                           method = "penalized_increasing")
      
      gamma_hat_lambda = result$par
      inc_gamma_hat_list[[index]] = gamma_hat_lambda
      gamma_exp_hat_lambda = c(gamma_hat_lambda[1:2],exp(gamma_hat_lambda[3:(m+1)]))
      
      theta_hat_lambda = S_mono %*% gamma_exp_hat_lambda
      theta_hat_lambda = as.numeric(theta_hat_lambda)
      names(theta_hat_lambda) = c("scale", rep("shape",m))
      inc_theta_hat_list[[index]] = theta_hat_lambda
      
      inc_neg_LL[index] = mono_neg_log_likelihood(gamma_hat_lambda, S_mono, x, q)        
      inc_penalized_neg_LL[index] = result$value
      
      xi_hat = theta_hat_lambda[-1]                                               
      sigma_hat = theta_hat_lambda[1] + cumsum(c(0, xi_hat[-m] * w[-m]))              
      phi_hat = xi_hat/sigma_hat
      psi_hat = c(sigma_hat[1],phi_hat)
      
      exp_info =(length(y) * exp.info.algebraic(psi_hat, x, q))
      svd = svd(exp_info + lambda * S)
      
      tryCatch({
        prod = solve(exp_info + lambda * S) %*% exp_info
        std_err = diag(solve(exp_info + lambda * S))
      }, error = function(e){
        U = svd$u
        L = diag(svd$d)
        V = svd$v
        L_inv = diag(1/svd$d)
        inv_svd = V%*%L_inv%*%t(U)
        prod = inv_svd%*%exp_info
        std_err = diag(inv_svd)
      })
      
      if(all(std_err >= 0)){
        std_err = sqrt(std_err)
      }else{
        std_err[which(std_err < 0)] = NA
        std_err[which(std_err >= 0)] = sqrt(std_err[which(std_err >= 0)])
      }
      
      inc_exp_info_list[[index]] = exp_info
      inc_svd_list[[index]] = svd
      inc_std_err_list[[index]] = std_err
      
      eigen_H = eigen(exp_info+lambda*S)
      
      inc_df[index] = sum(diag(prod))
      inc_AIC[index] = 2 * inc_neg_LL[index] + 2 * inc_df[index]
      inc_AICc[index] = inc_AIC[index]+(2 * inc_df[index] * (inc_df[index] + 1))/(n - inc_df[index] - 1)
      inc_BIC[index] = 2 * inc_neg_LL[index] + inc_df[index] * log(n)
      
      if(prod(svd$d)>.Machine$double.xmin){
        inc_marginal_neg_LL[index] = result$value + 0.5 * log(prod(svd$d)) - 0.5 * (m - 1) * log(lambda)
      }else if(det(exp_info+lambda*S)>.Machine$double.xmin){
        inc_marginal_neg_LL[index]=result$value + 0.5 * log(det(exp_info+lambda*S)) - 0.5 * (m - 1) * log(lambda)
      }else if(prod(eigen_H$values)>.Machine$double.xmin){
        inc_marginal_neg_LL[index]=result$value + 0.5 * log(prod(eigen_H$values)) - 0.5 * (m - 1) * log(lambda)
      }else{
        inc_marginal_neg_LL[index]=.Machine$double.xmax
      }
      
      gamma_loop = gamma_hat_lambda
    }
    inc_I_opt = min(which(-inc_marginal_neg_LL>=max(-inc_marginal_neg_LL)-0.01))
    #inc_I_opt = which.max(-inc_marginal_neg_LL)
    
    inc_theta_opt = inc_theta_hat_list[[inc_I_opt]] 
    inc_gamma_opt = inc_gamma_hat_list[[inc_I_opt]]
    
    lower = log10(val[max(inc_I_opt-2,1)])
    upper = log10(val[min(inc_I_opt+2,l)])
    
    inc_marg_neg_LL_optimized = optimize(marg_neg_LL,
                                         interval = c(lower,upper), 
                                         tol = .Machine$double.eps^0.25,
                                         maximum = FALSE,
                                         x = x, q = q, gamma_opt = inc_gamma_opt,
                                         method = "penalized_increasing")
    
    inc_lambda_best = 10^(inc_marg_neg_LL_optimized$minimum)
    
    result = optimize_neg_log_likelihood(mono_penalized_neg_log_likelihood, 
                                         gamma_opt = inc_gamma_opt,
                                         x = x, q = q, lambda = inc_lambda_best,
                                         method = "penalized_increasing")
    
    inc_gamma_best = result$par
    inc_gamma_exp_best = c(inc_gamma_best[1:2], exp(inc_gamma_best[3:(m+1)]))
    inc_theta_best = S_inc %*% inc_gamma_exp_best
    inc_theta_best = as.numeric(inc_theta_best)
    names(inc_theta_best) = c("scale", rep("shape",m))
    
    inc_neg_LL_best = mono_neg_log_likelihood(inc_gamma_best, S_inc, x, q)
    inc_penalized_neg_LL_best = result$value
    
    
    xi_hat = inc_theta_best[-1]                                                                # (xi_1,...,xi_m)
    sigma_hat = inc_theta_best[1] + cumsum(c(0, xi_hat[-m] * w[-m]))                                      # (sigma_1,...,sigma_m)
    phi_hat = xi_hat/sigma_hat
    psi_hat = c(sigma_hat[1], phi_hat)
    
    exp_info = length(y) * exp.info.algebraic(psi_hat, x, q)
    svd = svd(exp_info + inc_lambda_best * S)
    
    tryCatch({
      prod = solve(exp_info + inc_lambda_best * S) %*% exp_info
      std_err_best = diag(solve(exp_info + inc_lambda_best * S))
    },error = function(e){
      U = svd$u
      L = diag(svd$d)
      V = svd$v
      L_inv = diag(1/svd$d)
      inv_svd = V %*% L_inv %*% t(U)
      prod = inv_svd %*% exp_info
      std_err_best = diag(inv_svd)
    })
    
    if(all(std_err_best >= 0)){
      std_err_best = sqrt(std_err_best)
    }else{
      std_err_best[which(std_err_best < 0)] = NA
      std_err_best[which(std_err_best >= 0)] = sqrt(std_err_best[which(std_err_best >= 0)])
    }
    
    res$increasing$exp_info = exp_info
    res$increasing$svd_exp_info = svd
    
    tryCatch({
      res$increasing$std_err = std_err_best
    }, error=function(e){
      res$increasing$std_err = NA
    })
    
    eigen_H = eigen(exp_info + inc_lambda_best * S)

    
    inc_df_best = sum(diag(prod))
    inc_AIC_best = 2 * inc_neg_LL_best + 2 * inc_df_best
    inc_AICc_best= inc_AIC_best + (2 * inc_df_best * (inc_df_best + 1))/(n - inc_df_best - 1)
    inc_BIC_best = 2 * inc_neg_LL_best + inc_df_best * log(n)
    
    if(prod(svd$d)>.Machine$double.xmin){
      inc_marginal_neg_LL_best = result$value + 0.5 * log(prod(svd$d)) - 0.5 * (m - 1) * log(lambda)
    }else if(det(exp_info+inc_lambda_best*S)>.Machine$double.xmin){
      inc_marginal_neg_LL_best = result$value + 0.5 * log(det(exp_info + inc_lambda_best * S)) - 0.5 * (m - 1) * log(inc_lambda_best)
    }else if(prod(eigen_H$values)>.Machine$double.xmin){
      inc_marginal_neg_LL_best = result$value + 0.5 * log(prod(eigen_H$values)) - 0.5 * (m - 1) * log(inc_lambda_best)
    }else{
      inc_marginal_neg_LL_best = .Machine$double.xmax
    }
    
    # ---- Save optimal results ----
    
    res$increasing$theta = inc_theta_best
    res$increasing$gamma = inc_gamma_best
    res$increasing$neg_LL = inc_neg_LL_best
    res$increasing$penalized_neg_LL = inc_penalized_neg_LL_best
    res$increasing$marginal_neg_LL = inc_marginal_neg_LL_best
    res$increasing$AIC = inc_AIC_best
    res$increasing$AIC_c = inc_AICc_best
    res$increasing$BIC = inc_BIC_best
    res$increasing$lambda = inc_lambda_best
    res$increasing$df = inc_df_best
    
    # ---- Save other results ----
    
    res$increasing$AIC_vec = inc_AIC
    res$increasing$AICc_vec = inc_AICc
    res$increasing$BIC_vec = inc_BIC
    res$increasing$df_vec = inc_df
    res$increasing$marginal_neg_LL_vec = inc_marginal_neg_LL
    res$increasing$neg_LL_vec = inc_neg_LL
    res$increasing$penalized_neg_LL_vec = inc_penalized_neg_LL
    res$increasing$theta_hat_list = inc_theta_hat_list
    res$increasing$exp_info_list = inc_exp_info_list
    res$increasing$svd_exp_info_list = inc_svd_list
    res$increasing$std_err_list = inc_std_err_list
  }
  
  # ---- CASE 2: MONOTONICALLY DECREASING SHAPE PARAMETERS ----
  if(penalized_decreasing == TRUE){
    
    matrix_of_ones = matrix(1, nrow = m, ncol = m)
    
    S_dec = - matrix_of_ones * lower.tri(matrix_of_ones, diag = TRUE)
    S_dec[, 1] = rep(1,m)
    S_dec = rbind(c(rep(0, m)), S_dec)
    S_dec = cbind(c(1, rep(0, m)), S_dec)
    
    Ip = diag(c(0, 0, rep(1, m - 1)))
    
    # ---- Initialize vectors ----
    index = 0
    dec_penalized_neg_LL = rep(NA,l)
    dec_neg_LL = rep(NA,l)
    
    dec_df = rep(NA,l)
    dec_AIC = rep(NA,l)
    dec_AICc = rep(NA,l)
    dec_BIC = rep(NA,l)
    dec_marginal_neg_LL= rep(NA,l)
    
    dec_theta_hat_list = list()
    dec_gamma_hat_list = list()
    dec_exp_info_list = list()
    dec_svd_list = list()
    dec_std_err_list = list()
    
    sgn = -1  # sign - for decreasing (+ for increasing)
    S_mono = S_dec
    gamma_init = c(gamma_init[1], c(1.5, -c(1:(m-1))))
    gamma_loop = gamma_init
    
    # ---- "for" loop over lambda(s) ----
    for(lambda in val){
      
      index = index + 1
      
      result = optimize_neg_log_likelihood(mono_penalized_neg_log_likelihood,
                                           gamma_opt = gamma_loop,
                                           x = x, q = q, lambda = lambda,
                                           method  = "penalized_decreasing")
      
      gamma_hat_lambda = result$par
      dec_gamma_hat_list[[index]] = gamma_hat_lambda
      gamma_exp_hat_lambda = c(gamma_hat_lambda[1:2], exp(gamma_hat_lambda[3:(m+1)]))
      
      theta_hat_lambda = S_dec %*% gamma_exp_hat_lambda
      theta_hat_lambda = as.numeric(theta_hat_lambda)
      names(theta_hat_lambda) = c("scale", rep("shape", m))
      dec_theta_hat_list[[index]] = theta_hat_lambda
      
      dec_neg_LL[index] = mono_neg_log_likelihood(gamma_hat_lambda, S_dec, x, q)      
      dec_penalized_neg_LL[index] = result$value
      
      xi_hat = theta_hat_lambda[-1]                                                              
      sigma_hat = theta_hat_lambda[1] + cumsum(c(0,xi_hat[-m] * w[-m]))                                   
      phi_hat = xi_hat/sigma_hat
      psi_hat = c(sigma_hat[1], phi_hat)
      
      exp_info = (length(y) * exp.info.algebraic(psi_hat, x, q))
      svd = svd(exp_info + lambda * S)
      
      tryCatch({
        prod = solve(exp_info + lambda * S) %*% exp_info
        std_err = diag(solve(exp_info + lambda * S))
      }, error = function(e){
        U = svd$u
        L = diag(svd$d)
        V = svd$v
        L_inv = diag(1/svd$d)
        inv_svd = V %*% L_inv %*% t(U)
        prod = inv_svd %*% exp_info
        std_err = diag(inv_svd)
      })
      
      if(all(std_err >= 0)){
        std_err = sqrt(std_err)
      }else{
        std_err[which(std_err < 0)] = NA
        std_err[which(std_err >= 0)] = sqrt(std_err[which(std_err >= 0)])
      }
      
      dec_exp_info_list[[index]] = exp_info
      dec_svd_list[[index]] = svd
      dec_std_err_list[[index]] = std_err
      
      eigen_H = eigen(exp_info + lambda * S)
      
      dec_df[index] = sum(diag(prod))
      dec_AIC[index] = 2 * dec_neg_LL[index] + 2 * dec_df[index]
      dec_AICc[index] = dec_AIC[index] + (2 * dec_df[index] * (dec_df[index] + 1))/(n - dec_df[index] - 1)
      dec_BIC[index] = 2 * dec_neg_LL[index] + dec_df[index] * log(n)
      
      if(prod(svd$d)>.Machine$double.xmin){
        dec_marginal_neg_LL[index] = result$value + 0.5 * log(prod(svd$d)) - 0.5 * (m - 1) * log(lambda)
      }else if(det(exp_info+lambda*S)>.Machine$double.xmin){
        dec_marginal_neg_LL[index] = result$value + 0.5 * log(det(exp_info + lambda * S)) - 0.5 * (m - 1) * log(lambda)
      }else if(prod(eigen_H$values)>.Machine$double.xmin){
        dec_marginal_neg_LL[index] = result$value + 0.5 * log(prod(eigen_H$values)) - 0.5 * (m - 1) * log(lambda)
      }else{
        dec_marginal_neg_LL[index] = .Machine$double.xmax
      }
      
      gamma_loop = gamma_hat_lambda
    }
    
    # ---- Optimal penalties ----
    dec_I_opt = min(which(-dec_marginal_neg_LL >= max(-dec_marginal_neg_LL) - 0.01))
    #dec_I_opt = which.max(-dec_marginal_neg_LL)
    
    dec_theta_opt = dec_theta_hat_list[[dec_I_opt]] 
    dec_gamma_opt = dec_gamma_hat_list[[dec_I_opt]]
    
    
    lower = log10(val[max(dec_I_opt - 2, 1)])
    upper = log10(val[min(dec_I_opt + 2, l)])
    
    
    dec_marg_neg_LL_optimized = optimize(marg_neg_LL,
                                         interval = c(lower, upper), 
                                         tol = .Machine$double.eps^0.25,
                                         maximum = FALSE,
                                         x = x, q = q, gamma_opt = dec_gamma_opt,
                                         method = "penalized_decreasing")
    
    dec_lambda_best = 10^(dec_marg_neg_LL_optimized$minimum)
    
    result = optimize_neg_log_likelihood(mono_penalized_neg_log_likelihood, 
                                         gamma_opt = dec_gamma_opt,
                                         x = x, q = q, lambda = dec_lambda_best,
                                         method = "penalized_decreasing")
    
    dec_gamma_best = result$par
    dec_gamma_exp_best = c(dec_gamma_best[1:2], exp(dec_gamma_best[3:(m+1)]))
    dec_theta_best = S_dec %*% dec_gamma_exp_best
    dec_theta_best = as.numeric(dec_theta_best)
    names(dec_theta_best) = c("scale", rep("shape", m))
    
    dec_neg_LL_best = mono_neg_log_likelihood(dec_gamma_best, S_dec, x, q)
    dec_penalized_neg_LL_best = result$value
    
    xi_hat = dec_theta_best[-1]                                                                # (xi_1,...,xi_m)
    sigma_hat = dec_theta_best[1] + cumsum(c(0, xi_hat[-m] * w[-m]))                                      # (sigma_1,...,sigma_m)
    phi_hat = xi_hat/sigma_hat
    psi_hat = c(sigma_hat[1], phi_hat)
    
    exp_info = length(y) * exp.info.algebraic(psi_hat, x, q)
    svd = svd(exp_info + dec_lambda_best * S)
    
    tryCatch({
      prod = solve(exp_info+dec_lambda_best * S) %*% exp_info
      std_err_best = diag(solve(exp_info + dec_lambda_best * S))
    },error = function(e){
      U = svd$u
      L = diag(svd$d)
      V = svd$v
      L_inv = diag(1/svd$d)
      inv_svd = V %*% L_inv %*% t(U)
      prod = inv_svd %*% exp_info
      std_err_best = diag(inv_svd)
    })
    
    if(all(std_err_best >= 0)){
      std_err_best = sqrt(std_err_best)
    }else{
      std_err_best[which(std_err_best < 0)] = NA
      std_err_best[which(std_err_best >= 0)] = sqrt(std_err_best[which(std_err_best >= 0)])
    }
    
    res$decreasing$exp_info = exp_info
    res$decreasing$svd_exp_info = svd
    
    tryCatch({
      res$decreasing$std_err = std_err_best
    }, error=function(e){
      res$decreasing$std_err = NA
    })
    
    #res$decreasing$std_err = std_err_best
    
    eigen_H = eigen(exp_info + dec_lambda_best * S)
    
    dec_df_best = sum(diag(prod))
    dec_AIC_best = 2 * dec_neg_LL_best + 2 * dec_df_best
    dec_AICc_best = dec_AIC_best + (2 * dec_df_best * (dec_df_best + 1))/(n - dec_df_best - 1)
    dec_BIC_best = 2 * dec_neg_LL_best + dec_df_best * log(n)
    
    if(prod(svd$d)>.Machine$double.xmin){
      dec_marginal_neg_LL_best = result$value + 0.5 * log(prod(svd$d)) - 0.5 * (m - 1) * log(lambda)
    }else if(det(exp_info+dec_lambda_best*S)>.Machine$double.xmin){
      dec_marginal_neg_LL_best = result$value + 0.5 * log(det(exp_info + dec_lambda_best * S)) - 0.5 * (m - 1) * log(dec_lambda_best)
    }else if(prod(eigen_H$values)>.Machine$double.xmin){
      dec_marginal_neg_LL_best = result$value + 0.5 * log(prod(eigen_H$values)) - 0.5 * (m - 1) * log(dec_lambda_best)
    }else{
      dec_marginal_neg_LL_best = .Machine$double.xmax
    }
    
    # ---- Save optimal results ----
    
    res$decreasing$theta = dec_theta_best
    res$decreasing$gamma = dec_gamma_best
    res$decreasing$neg_LL = dec_neg_LL_best
    res$decreasing$penalized_neg_LL = dec_penalized_neg_LL_best
    res$decreasing$marginal_neg_LL = dec_marginal_neg_LL_best
    res$decreasing$AIC = dec_AIC_best
    res$decreasing$AICc = dec_AICc_best
    res$decreasing$BIC = dec_BIC_best
    res$decreasing$lambda = dec_lambda_best
    res$decreasing$df = dec_df_best
    
    # ---- Save other results ----
    
    res$decreasing$AIC_vec = dec_AIC
    res$decreasing$AICc_vec = dec_AICc
    res$decreasing$BIC_vec = dec_BIC
    res$decreasing$df_vec = dec_df
    res$decreasing$marginal_neg_LL_vec = dec_marginal_neg_LL
    res$decreasing$neg_LL_vec = dec_neg_LL
    res$decreasing$penalized_neg_LL_vec = dec_penalized_neg_LL
    res$decreasing$theta_hat_list = dec_theta_hat_list
    res$decreasing$exp_info_list = dec_exp_info_list
    res$decreasing$svd_exp_info_list = dec_svd_list
    res$decreasing$std_err_list = dec_std_err_list
  }
  
  return(res)
}
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
replicate_m.t.GP.fit = function(replications, n, q, params, distribution,
                                unpenalized = TRUE, penalized = TRUE,
                                penalized_weighted = FALSE,
                                penalized_increasing = FALSE, 
                                penalized_decreasing = FALSE,
                                val = c(10^seq(-2, 16, by = 1)),
                                save_result = FALSE){
  # function that fits a multiple-threshold GP model to 
  # multiple datasets generated from the same distribution
  
  # possible distributions are:
  # "exponential", "gp", "weibull", "half-normal", "half-cauchy", 
  # "hybrid uniform-gp", "hybrid uniform-weibull"
  
  B = replications
  num_cores = 6
  cl = makeCluster(num_cores)
  registerDoParallel(cl)
  getDoParWorkers()
  
  res = list()
  final_result = vector("list",B)
  
  writeLines(c(""), "log.txt")
  
  if(distribution == "exponential"){
    start_time = Sys.time()
    sim_result = foreach(b = 1:B) %dopar%{
      sink("log.txt", append=TRUE)
      cat(paste("Starting iteration",b, "Time elapsed", Sys.time()-start_time, "\n"))
      sink()
      tryCatch({
        set.seed(b)
        x = rexp(n, params)
        fit_true = MASS::fitdistr(x, "exponential")
        result = m.t.GP.fit(x, q,
                            unpenalized = unpenalized, 
                            penalized = penalized,
                            penalized_weighted = penalized_weighted,
                            penalized_increasing = penalized_increasing,
                            penalized_decreasing = penalized_decreasing, 
                            val = val)
        
        result$true_fit$estimates = fit_true$estimate
        result$true_fit$std_err = fit_true$sd
        result$true_fit$neg_LL = -fit_true$loglik
        final_result[[b]] = result
      }, error = function(e){
        final_result[[b]] = NA
      })
    }
    end_time = Sys.time()
    res = sim_result
    res$running_time = end_time - start_time
    
    stopCluster(cl)
    
    if(save_result == TRUE){
      title = paste0("EXP_rate_", params, "_replications_", B, ".rds")
      saveRDS(sim_result, title)
    }
    return(res)
  }
  
  if(distribution == "weibull"){
    start_time = Sys.time()
    sim_result = foreach(b = 1:B) %dopar%{
      sink("log.txt", append=TRUE)
      cat(paste("Starting iteration",b, "Time elapsed", Sys.time()-start_time, "\n"))
      sink()
      tryCatch({
        set.seed(b)
        x = rweibull(n, shape = params[1], scale = params[2])
        fit_true = MASS::fitdistr(x, "weibull")
        result = m.t.GP.fit(x, q, 
                            unpenalized = unpenalized, 
                            penalized = penalized,
                            penalized_weighted = penalized_weighted,
                            penalized_increasing = penalized_increasing,
                            penalized_decreasing = penalized_decreasing, 
                            val = val)
        
        result$true_fit$estimates = fit_true$estimate
        result$true_fit$std_err = fit_true$sd
        result$true_fit$neg_LL = -fit_true$loglik
        final_result[[b]] = result
      }, error = function(e){
        final_result[[b]] = NA
      })
    }
    end_time = Sys.time()
    res = sim_result
    res$running_time = end_time - start_time
    
    stopCluster(cl)
    
    if(save_result == TRUE){
      title = paste0("WEIBULL_shape_", params[1],
                        "_scale_", params[2], "_replications_", B, ".rds")
      saveRDS(sim_result, title)
    }
    return(res)
  }
  
  if(distribution == "half-normal"){
    start_time = Sys.time()
    sim_result = foreach(b = 1:B) %dopar%{
      sink("log.txt", append=TRUE)
      cat(paste("Starting iteration", b, "Time elapsed", Sys.time()-start_time, "\n"))
      sink()
      tryCatch({
        set.seed(b)
        x = rnorm(n, mean = params[1], sd = params[2])
        fit_true = MASS::fitdistr(x, "normal")
        x = abs(x)
        result = m.t.GP.fit(x, q, 
                            unpenalized = unpenalized, 
                            penalized = penalized,
                            penalized_weighted = penalized_weighted,
                            penalized_increasing = penalized_increasing,
                            penalized_decreasing = penalized_decreasing, 
                            val = val)
        
        result$true_fit$estimates = fit_true$estimate
        result$true_fit$std_err = fit_true$sd
        result$true_fit$neg_LL = -fit_true$loglik
        final_result[[b]] = result
      }, error = function(e){
        final_result[[b]] = NA
      })
    }
    end_time = Sys.time()
    res = sim_result
    res$running_time = end_time - start_time
    
    stopCluster(cl)
    
    if(save_result == TRUE){
      title = paste0("HalfNormal_mean_", params[1],
                     "_sd_", params[2], "_replications_", B, ".rds")
      saveRDS(sim_result, title)
    }
    return(res)
  }
  
  if(distribution == "half-cauchy"){
    start_time = Sys.time()
    sim_result = foreach(b = 1:B) %dopar%{
      sink("log.txt", append=TRUE)
      cat(paste("Starting iteration", b, "Time elapsed", Sys.time()-start_time, "\n"))
      sink()
      tryCatch({
        set.seed(b)
        x = rcauchy(n, location = params[1], scale = params[2])
        fit_true = MASS::fitdistr(x, "cauchy")
        x = abs(x)
        result = m.t.GP.fit(x, q, 
                            unpenalized = unpenalized, 
                            penalized = penalized,
                            penalized_weighted = penalized_weighted,
                            penalized_increasing = penalized_increasing,
                            penalized_decreasing = penalized_decreasing, 
                            val = val)
        
        result$true_fit$estimates = fit_true$estimate
        result$true_fit$std_err = fit_true$sd
        result$true_fit$neg_LL = -fit_true$loglik
        final_result[[b]] = result
      }, error = function(e){
        final_result[[b]] = NA
      })
    }
    end_time = Sys.time()
    res = sim_result
    res$running_time = end_time - start_time
    
    stopCluster(cl)
    
    if(save_result == TRUE){
      title = paste0("HalfCauchy_location_", params[1],
                     "_scale_", params[2], "_replications_", B, ".rds")
      saveRDS(sim_result, title)
    }
    return(res)
  }
 
  if(distribution == "gp"){
    start_time = Sys.time()
    sim_result = foreach(b = 1:B) %dopar%{
      sink("log.txt", append=TRUE)
      cat(paste("Starting iteration", b, "Time elapsed", Sys.time()-start_time, "\n"))
      sink()
      tryCatch({
        set.seed(b)
        x = mev::rgp(n, scale = params[1], shape = params[2]) 
        fit_true = mev::fit.gpd(x)
        result = m.t.GP.fit(x, q,
                            unpenalized = unpenalized, 
                            penalized = penalized,
                            penalized_weighted = penalized_weighted,
                            penalized_increasing = penalized_increasing,
                            penalized_decreasing = penalized_decreasing, 
                            val = val)
        
        result$true_fit$estimates = fit_true$estimate
        result$true_fit$std_err = fit_true$std.err
        result$true_fit$neg_LL = fit_true$nllh
        
        final_result[[b]] = result
      }, error = function(e){
        final_result[[b]] = NA
      })
    }
    end_time = Sys.time()
    res = sim_result
    res$running_time = end_time - start_time
    
    stopCluster(cl)
    
    if(save_result == TRUE){
      title = paste0("GP_scale_", params[1], "_shape_", params[2], 
                     "_replications_", B, ".rds")
      saveRDS(sim_result, title)
    }
    return(res)
  }
  
  if(distribution == "hybrid uniform-gp"){
    start_time = Sys.time()
    sim_result = foreach(b = 1:B) %dopar%{
      sink("log.txt", append=TRUE)
      cat(paste("Starting iteration", b, "Time elapsed", Sys.time()-start_time, "\n"))
      sink()
      tryCatch({
        set.seed(b)
        x = r_hybrid_unif_gp(n, params[1], params[2], params[3], params[4], params[5])
        fit_true = fit_hybrid_unif_gp(x)
        result = m.t.GP.fit(x, q,
                            unpenalized = unpenalized, 
                            penalized = penalized,
                            penalized_weighted = penalized_weighted,
                            penalized_increasing = penalized_increasing,
                            penalized_decreasing = penalized_decreasing, 
                            val = val)
        
        result$true_fit$estimates = fit_true$estimate
        result$true_fit$std_err = fit_true$sd
        result$true_fit$neg_LL = -fit_true$loglik
        
        final_result[[b]] = result
      }, error = function(e){
        final_result[[b]] = NA
      })
    }
    end_time = Sys.time()
    res = sim_result
    res$running_time = end_time - start_time
    
    stopCluster(cl)
    
    if(save_result == TRUE){
      title = paste0("Hybrid_Unif_GP_lower_", params[1], "_upper_", params[2],
                     "_shape_", params[3], "_scale_", params[4], 
                     "_mixing_parameter_", params[5],
                     "_replications_", B, ".rds")
      saveRDS(sim_result, title)
    }
    return(res)
  }
  
  if(distribution == "hybrid uniform-weibull"){
    start_time = Sys.time()
    sim_result = foreach(b = 1:B) %dopar%{
      sink("log.txt", append=TRUE)
      cat(paste("Starting iteration", b, "Time elapsed", Sys.time()-start_time, "\n"))
      sink()
      tryCatch({
        set.seed(b)
        x = r_hybrid_unif_weibull(n, params[1], params[2], params[3], params[4], params[5])
        fit_true = fit_hybrid_unif_weibull(x)
        result = m.t.GP.fit(x, q,
                            unpenalized = unpenalized, 
                            penalized = penalized,
                            penalized_weighted = penalized_weighted,
                            penalized_increasing = penalized_increasing,
                            penalized_decreasing = penalized_decreasing, 
                            val = val)
        
        result$true_fit$estimates = fit_true$estimate
        result$true_fit$std_err = fit_true$sd
        result$true_fit$neg_LL = -fit_true$loglik
        
        final_result[[b]] = result
      }, error = function(e){
        final_result[[b]] = NA
      })
    }
    end_time = Sys.time()
    res = sim_result
    res$running_time = end_time - start_time
    
    stopCluster(cl)
    
    if(save_result == TRUE){
      title = paste0("Hybrid_Unif_Weibull_lower_", params[1], "_upper_", params[2],
                     "_shape_", params[3], "_scale_", params[4], 
                     "_mixing_parameter_", params[5],
                     "_replications_", B, ".rds")
      saveRDS(sim_result, title)
    }
    return(res)
  }
}
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
compute_quantiles_m.t.GP.fit = function(alpha, data, q, params, 
                                        distribution, repl,
                                        unpenalized = TRUE, penalized = TRUE,
                                        penalized_weighted = FALSE,
                                        penalized_increasing = FALSE, 
                                        penalized_decreasing = FALSE,
                                        val = c(10^seq(-2, 16, by = 1)),
                                        save_result = FALSE){
  
  # function that computes quantiles at a vector of probabilities alpha, given 
  # the data (and corresponding estimates from the m.t. GP fits), and parameters 
  # of the true underlying distribution that generated the data
  
  # quantiles are computed for several replications of the data
  
  l = length(alpha)
  
  q_inc_pen_gpd = rep(NA,l)
  q_dec_pen_gpd = rep(NA,l)
  q_w_pen_gpd = rep(NA,l)
  q_pen_gpd = rep(NA,l)
  q_unpen_gpd = rep(NA,l)
  q_true_fit = rep(NA,l)
  q_true = rep(NA,l)
  
  B = repl
  
  num_cores = 6
  cl = makeCluster(num_cores)
  registerDoParallel(cl)
  getDoParWorkers()
  
  final_result = vector("list", B)
  
  writeLines(c(""), "log.txt")
  
  
  start_time = Sys.time()
  quant_result = foreach(b = 1:B, .export=c("qgpd_multiple_thresholds")) %dopar% {
    sink("log.txt", append=TRUE)
    cat(paste("Starting iteration",b, "Time elapsed", Sys.time()-start_time, "\n"))
    sink()
    result = list()
    #tryCatch({
      x = data[[b]]$info$data
      
      if(unpenalized == TRUE){
        for(i in 1:l){
          q_unpen_gpd[i] = qgpd_multiple_thresholds(x, q, alpha[i], data[[b]]$unpenalized$theta)
        }
        result$q_unpen_gpd = q_unpen_gpd
      }
      
      if(penalized == TRUE){
        for(i in 1:l){
          q_pen_gpd[i] = qgpd_multiple_thresholds(x, q, alpha[i], data[[b]]$penalized$theta)
        }
        result$q_pen_gpd = q_pen_gpd
      }

      if(penalized_weighted == TRUE){
        for(i in 1:l){
          q_w_pen_gpd[i] = qgpd_multiple_thresholds(x, q, alpha[i], data[[b]]$weighted$theta)
        }
        result$q_w_pen_gpd = q_w_pen_gpd
      }

      if(penalized_increasing == TRUE){
        for(i in 1:l){
          q_inc_pen_gpd[i] = qgpd_multiple_thresholds(x, q, alpha[i], data[[b]]$increasing$theta)
        }
        result$q_inc_pen_gpd = q_inc_pen_gpd
      }

      if(penalized_decreasing == TRUE){
        for(i in 1:l){
          q_dec_pen_gpd[i] = qgpd_multiple_thresholds(x, q, alpha[i], data[[b]]$decreasing$theta)
        }
        result$q_dec_pen_gpd = q_dec_pen_gpd
      }

      if(distribution == "exponential"){
        for(i in 1:l){
          q_true_fit[i] = qexp(alpha[i], data[[b]]$true_fit$estimates)
          q_true[i] = qexp(alpha[i], params)
        }
        result$q_true_fit = q_true_fit
        result$q_true = q_true
      }

      if(distribution == "weibull"){
        for(i in 1:l){
          q_true_fit[i] = qweibull(alpha[i], shape = data[[b]]$true_fit$estimates[1],
                                   scale = data[[b]]$true_fit$estimates[2])
          q_true[i] = qweibull(alpha[i], shape = params[1], scale = params[2])
        }
        result$q_true_fit = q_true_fit
        result$q_true = q_true
      }

      if(distribution == "half-normal"){
        for(i in 1:l){
          q_true_fit[i] = qnorm((1 + alpha[i])/2,
                                mean = data[[b]]$true_fit$estimates[1],
                                sd = data[[b]]$true_fit$estimates[2])
          q_true[i] = qnorm((1+alpha[i])/2,mean=mu,sd=sig)
        }
        result$q_true_fit = q_true_fit
        result$q_true = q_true
      }

      if(distribution == "half-cauchy"){
        for(i in 1:l){
          q_true_fit[i] = qcauchy((1 + alpha[i])/2,
                                  location = data[[b]]$true_fit$estimates[1],
                                  scale = data[[b]]$true_fit$estimates[2])
          q_true[i] = qcauchy((1 + alpha[i])/2, location = 0, scale = 1)
        }
        result$q_true_fit = q_true_fit
        result$q_true = q_true
      }

      if(distribution == "gp"){
        for(i in 1:l){
          q_true_fit[i] = mev::qgp(alpha[i],
                                   scale = data[[b]]$true_fit$estimates[1],
                                   shape = data[[b]]$true_fit$estimates[2])
          q_true[i] = mev::qgp(alpha[i], scale = params[1], shape = params[2])
        }
        result$q_true_fit = q_true_fit
        result$q_true = q_true
      }

      if(distribution == "hybrid uniform-gp"){
        for(i in 1:l){
          q_true_fit[i] = q_hybrid_unif_gp(alpha[i],
                                           data[[b]]$true_fit$estimates[1],
                                           data[[b]]$true_fit$estimates[2],
                                           data[[b]]$true_fit$estimates[3],
                                           data[[b]]$true_fit$estimates[4],
                                           data[[b]]$true_fit$estimates[5])
          q_true[i] = q_hybrid_unif_gp(alpha[i], params[1], params[2],
                                       params[3], params[4], params[5])
        }
        result$q_true_fit = q_true_fit
        result$q_true = q_true
      }

      if(distribution == "hybrid uniform-weibull"){
        for(i in 1:l){
          q_true_fit[i] = q_hybrid_unif_weibull(alpha[i],
                                                data[[b]]$true_fit$estimates[1],
                                                data[[b]]$true_fit$estimates[2],
                                                data[[b]]$true_fit$estimates[3],
                                                data[[b]]$true_fit$estimates[4],
                                                data[[b]]$true_fit$estimates[5])
          q_true[i] = q_hybrid_unif_weibull(alpha[i], params[1], params[2],
                                            params[3], params[4], params[5])
        }
        result$q_true_fit = q_true_fit
        result$q_true = q_true
      }
      
      final_result[[b]] = result
      
    #}, error = function(e){
    #  final_result[[b]] = NA
    #})
  }
  end_time = Sys.time()
  #res = quant_result
  #res$running_time = end_time - start_time   
  stopCluster(cl)
  
  if(save_result == TRUE){
    if(distribution == "exponential"){
      title = paste0("quantiles_EXP_rate_", params, "_replications_", B, ".rds")
    }

    if(distribution == "weibull"){
      title = paste0("quantiles_Weibull_shape_", params[1],
                     "_scale_", params[2], "_replications_", B, ".rds")
    }

    if(distribution == "half-normal"){
      title = paste0("quantiles_HalfNormal_mean_", params[1],
                     "_sd_", params[2], "_replications_", B, ".rds")
    }

    if(distribution == "half-cauchy"){
      title = paste0("quantiles_HalfCauchy_location_", params[1],
                     "_scale_", params[2], "_replications_", B, ".rds")
    }

    if(distribution == "gp"){
      title = paste0("quantiles_GP_scale_", params[1],
                     "_shape_", params[2], "_replications_", B, ".rds")
    }

    if(distribution == "hybrid uniform-gp"){
      title = paste0("quantiles_Uniform_GP_lower_", params[1],
                     "_upper_", params[2], "_shape_", params[3],
                     "_scale_", params[4], "_mixing_parameter_", params[5],
                     "_replications_", B, ".rds")
    }

    if(distribution == "hybrid uniform-weibull"){
      title = paste0("quantiles_Uniform_Weibull_lower_", params[1],
                     "_upper_", params[2], "_shape_", params[3],
                     "_scale_", params[4], "_mixing_parameter_", params[5],
                     "_replications_", B, ".rds")
    }

    saveRDS(quant_result, title)
  }
  
  return(quant_result)
}
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
compute_ise_mse_var_bias = function(quant_result, alpha, replications,
                                unpenalized = TRUE, penalized = TRUE,
                                penalized_weighted = FALSE,
                                penalized_increasing = FALSE, 
                                penalized_decreasing = FALSE,
                                val = c(10^seq(-2, 16, by = 1)),
                                save_result = FALSE){
  
  # function that computes ISE, MSE, VAR and BIAS of quantiles evaluated at a 
  # vector of probabilities alpha, given the quantile estimates obtained
  # from the m.t. GP fits
  
  # quantile estimates are given for several replications of the data
  
  l = length(alpha)
  B = replications
  
  res = list()
  
  q_true_fit = matrix(NA, nrow=B, ncol=l)
  q_true = matrix(NA, nrow=B, ncol=l)
  
  for(b in 1:B){
    if(!is.na(quant_result[[b]])[1]){
      q_true_fit[b, ] = quant_result[[b]]$q_true_fit
      q_true[b, ] = quant_result[[b]]$q_true
    }
  }
  q_true_fit = q_true_fit[complete.cases(q_true_fit),]
  q_true = q_true[complete.cases(q_true),]
  MSE_true_fit = colMeans((q_true - q_true_fit)^2)
  ISE_true_fit = rowMeans((q_true_fit - q_true)^2)
  VAR_true_fit = colMeans((t(t(q_true_fit)-colMeans(q_true_fit)))^2)
  BIAS_true_fit = colMeans(q_true_fit)-q_true[1,]
  res$MSE_true_fit = MSE_true_fit
  res$ISE_true_fit = ISE_true_fit
  res$VAR_true_fit = VAR_true_fit
  res$BIAS_true_fit = BIAS_true_fit
  
  
  if(unpenalized == TRUE){
    q_unpen_gpd = matrix(NA, nrow=B, ncol=l)
    for(b in 1:B){
      if(!is.na(quant_result[[b]])[1]){
        q_unpen_gpd[b,] = quant_result[[b]]$q_unpen_gpd
      }
    }
    q_unpen_gpd = q_unpen_gpd[complete.cases(q_unpen_gpd),]
    MSE_unpen_gpd = colMeans((q_unpen_gpd - q_true)^2)
    ISE_unpen_gpd = rowMeans((q_unpen_gpd - q_true)^2)
    VAR_unpen_gpd = colMeans((t(t(q_unpen_gpd) - colMeans(q_unpen_gpd)))^2)
    BIAS_unpen_gpd = colMeans(q_unpen_gpd) - q_true[1,]
    res$MSE_unpen_gpd = MSE_unpen_gpd
    res$ISE_unpen_gpd = ISE_unpen_gpd
    res$VAR_unpen_gpd = VAR_unpen_gpd
    res$BIAS_unpen_gpd = BIAS_unpen_gpd
  }
  
  if(penalized == TRUE){
    q_pen_gpd = matrix(NA, nrow=B, ncol=l)
    for(b in 1:B){
      if(!is.na(quant_result[[b]])[1]){
        q_pen_gpd[b,] = quant_result[[b]]$q_pen_gpd
      }
    }
    q_pen_gpd = q_pen_gpd[complete.cases(q_pen_gpd),]
    MSE_pen_gpd = colMeans((q_pen_gpd - q_true)^2)
    ISE_pen_gpd = rowMeans((q_pen_gpd - q_true)^2)
    VAR_pen_gpd = colMeans((t(t(q_pen_gpd) - colMeans(q_pen_gpd)))^2)
    BIAS_pen_gpd = colMeans(q_pen_gpd) - q_true[1,]
    res$MSE_pen_gpd = MSE_pen_gpd
    res$ISE_pen_gpd = ISE_pen_gpd
    res$VAR_pen_gpd = VAR_pen_gpd
    res$BIAS_pen_gpd = BIAS_pen_gpd
  }
  
  if(penalized_weighted == TRUE){
    q_w_pen_gpd = matrix(NA, nrow=B, ncol=l)
    for(b in 1:B){
      if(!is.na(quant_result[[b]])[1]){
        q_w_pen_gpd[b,] = quant_result[[b]]$q_w_pen_gpd
      }
    }
    q_w_pen_gpd = q_w_pen_gpd[complete.cases(q_w_pen_gpd),]
    MSE_w_pen_gpd = colMeans((q_w_pen_gpd - q_true)^2)
    ISE_w_pen_gpd = rowMeans((q_w_pen_gpd - q_true)^2)
    VAR_w_pen_gpd = colMeans((t(t(q_w_pen_gpd) - colMeans(q_w_pen_gpd)))^2)
    BIAS_w_pen_gpd = colMeans(q_w_pen_gpd) - q_true[1,]
    res$MSE_w_pen_gpd = MSE_w_pen_gpd
    res$ISE_w_pen_gpd = ISE_w_pen_gpd
    res$VAR_w_pen_gpd = VAR_w_pen_gpd
    res$BIAS_w_pen_gpd = BIAS_w_pen_gpd
  }
  
  if(penalized_increasing == TRUE){
    q_inc_pen_gpd = matrix(NA, nrow=B, ncol=l)
    for(b in 1:B){
      if(!is.na(quant_result[[b]])[1]){
        q_inc_pen_gpd[b,] = quant_result[[b]]$q_inc_pen_gpd
      }
    }
    q_inc_pen_gpd = q_inc_pen_gpd[complete.cases(q_inc_pen_gpd),]
    MSE_inc_pen_gpd = colMeans((q_inc_pen_gpd - q_true)^2)
    ISE_inc_pen_gpd = rowMeans((q_inc_pen_gpd - q_true)^2)
    VAR_inc_pen_gpd = colMeans((t(t(q_inc_pen_gpd) - colMeans(q_inc_pen_gpd)))^2)
    BIAS_inc_pen_gpd = colMeans(q_inc_pen_gpd) - q_true[1,]
    res$MSE_inc_pen_gpd = MSE_inc_pen_gpd
    res$ISE_inc_pen_gpd = ISE_inc_pen_gpd
    res$VAR_inc_pen_gpd = VAR_inc_pen_gpd
    res$BIAS_inc_pen_gpd = BIAS_inc_pen_gpd
  }
  
  
  if(penalized_decreasing == TRUE){
    q_dec_pen_gpd = matrix(NA, nrow=B, ncol=l)
    for(b in 1:B){
      if(!is.na(quant_result[[b]])[1]){
        q_dec_pen_gpd[b,] = quant_result[[b]]$q_dec_pen_gpd
      }
    }
    q_dec_pen_gpd = q_dec_pen_gpd[complete.cases(q_dec_pen_gpd),]
    MSE_dec_pen_gpd = colMeans((q_dec_pen_gpd - q_true)^2)
    ISE_dec_pen_gpd = rowMeans((q_dec_pen_gpd - q_true)^2)
    VAR_dec_pen_gpd = colMeans((t(t(q_dec_pen_gpd) - colMeans(q_dec_pen_gpd)))^2)
    BIAS_dec_pen_gpd = colMeans(q_dec_pen_gpd) - q_true[1,]
    res$MSE_dec_pen_gpd = MSE_dec_pen_gpd
    res$ISE_dec_pen_gpd = ISE_dec_pen_gpd
    res$VAR_dec_pen_gpd = VAR_dec_pen_gpd
    res$BIAS_dec_pen_gpd = BIAS_dec_pen_gpd
  }
  
  return(res)

}
# ------------------------------------------------------------------------------

# ==============================================================================


