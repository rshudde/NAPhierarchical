### This is the file with the forms for MLl of the Bayes factor testing
rm(list = ls())
library(invgamma)
library(mvtnorm)
library(BayesFactor)
library(testit)
library(plotrix)
library(magic)
library(doParallel)
# source("./RejectionSampler.R")
# source('./JZS.r')
# blah

## to help visuMLize data
rotate = function(x)
{
  rev = t(apply(x, 2, rev)) 
  return(rev)
}

## for cMLculating gamma(n) where n is large
sterling_gamma = function(n)
{
  const = n*log(n) - n
  return(const)
}

get_g = function(p, w)
{
  to_return = w/(p + 2)
  return(to_return)
}

###### functions for data generation 
get_sigma2 = function(alpha = 0.1, lambda = 0.1)
{
  sigma2 = rinvgamma(1, alpha, lambda)
  return(sigma2)
}

get_Cj = function(Xj)
{
  Cj = list()
  for (i in 1:length(Xj))
  {
    temp = Xj[[i]]
    Cj[[i]] = round(solve(crossprod(temp)),5)
  }
  return(Cj)
}

# J = number of groups
# p = number of covariates in each group
get_data = function(J, rho, p = 2, prop_mean = 0, sample_min = 10, sample_max = 10, tau2 = 0.09/2, generate_local = FALSE)
{
  # check t oavoid singularity for solve(X^TX)
  if (sample_min <= p) warning("make sure minimum n <= p")
  
  # opening paramters
  if (sample_min == sample_max) {
    nj = rep(sample_min, J)
  } else {
    nj = sample(sample_min:sample_max, J, replace = T)
  }
  
  B = diag(p)
  Cj = list()
  Xj = list()
  beta = list()
  a = 0 # mean for alpha draws
  
  # generate Xj matrix and Cj matrix, store as a list of matrices
  for (i in 1:J)
  {
    temp = matrix(rnorm(nj[i]*p), nrow = nj[i], ncol = p)
    temp = scale(temp)
    temp = sqrt(nj[i]/(nj[i] - 1))*temp
    Cj[[i]] = round(solve(crossprod(temp)),5)
    Xj[[i]] = temp
  }
  
  # generate parameters sigma2 and alpha
  # sigma2 = get_sigma2(1, 1)
  sigma2 = 1
  # alpha = mvtnorm::rmvnorm(n = 1, mean = rep(a, p), sigma = sigma2*B)
  alpha = numeric(p)
  
  # generate beta from local and non-local case
  beta = list()
  # if (generate_local)
  # {
  #   for (i in 1:J)
  #   {
  #     temp_sigma = tau2*sigma2*nj[i]*Cj[[i]]
  #     temp = c(rmvn(1, mu = alpha, Sigma = temp_sigma))
  #     beta[[i]] = temp
  #   }
  # } else {
  #   for (i in 1:J)
  #   {
  #     temp = rejection_sampler(Xj[[i]], tau2, sigma2 = sigma2, prop_mean = 1)$beta
  #     beta[[i]] = temp
  #   }
  # }
  
  #  draw correlated betas
  sigmaBeta2 = .001
  for (i in 1:J)
  {
    corr = matrix(rho, nrow = p, ncol = p)
    diag(corr) = 1
    # temp = c(t(mvtnorm::rmvnorm(n = 1, mean = c(-prop_mean, prop_mean), sigma = (.1^2)*corr)))
    temp = c(t(mvtnorm::rmvnorm(n = 1, mean = c(-prop_mean, prop_mean), sigma = sigmaBeta2*corr)))
    
    beta[[i]] = temp
  }
  
  
  # generate y values 
  yj = list()
  for (i in 1:J)
  {
    # temp = rmvn(1, mu = c(Xj[[i]] %*% beta[[i]]), Sigma = diag(sigma2, nrow = nj[i], ncol = nj[i]))
    temp = rnorm(nj[i], c(Xj[[i]] %*% beta[[i]]), sqrt(sigma2))
    yj[[i]] = c(temp)
  }
  
  return(list(yj = yj,
              betaj = beta,
              sigma2 = sigma2,
              alpha = alpha,
              Xj = Xj))
}

# null marginal - not sure it's of any use 
get_m0_OLD = function(data, tau2 = 0.09/2)
{
  Xj = data$Xj
  yj = data$yj
  betaj = data$betaj
  nj = unlist(lapply(Xj, nrow))
  N = sum(nj)
  J = length(Xj)
  p = ncol(Xj[[1]])
  B = diag(p)
  Binv = solve(B)
  Cj = get_Cj(Xj)
  a = rep(0, p)
  
  # starting definitinos for G and g
  yty = lapply(yj, crossprod)
  XtX = lapply(Xj, crossprod)
  Xty = matrix(0, nrow = ncol(Xj[[1]]), ncol = 1)
  for (i in 1:J)
  {
    Xty = Xty + crossprod(Xj[[i]], yj[[i]])
  }
  
  XtX = Reduce('+', XtX) # sum of Xj^TXj
  yty = sum(unlist(yty))
  
  # quantities needed
  Q = XtX + solve(B)
  b = Xty + solve(B) %*% a
  M = yty + crossprod(a, solve(B)) %*% a - crossprod(b, solve(Q)) %*% b
  
  const = (2*pi)^(N/2)*(det(B)*det(Q))^(-1/2)*sterling_gamma(N/2)
  quad = ((1/2)*M)^(-N/2)
  
  marginal = log(const*quad)
  to_return = as.numeric(marginal)
    
  return(to_return)
}

get_mL = function(data, tau2 = 0.09/2, a = rep(0, p), B = diag(p), s = 3, r = 1)
{
  Xj = data$Xj
  yj = data$yj
  # terms needed
  nj = unlist(lapply(Xj, nrow))
  N = sum(nj)
  p = ncol(Xj[[1]])
  J = length(Xj)
  Cj = get_Cj(Xj)
  Binv = solve(B)
  
  # calculate D
  Dj = list()
  for (i in 1:length(Xj))
  {
    Dj[[i]] = diag(nj[[i]]) + tau2 * nj[[i]] * (Xj[[i]] %*% Cj[[i]] %*% t(Xj[[i]]))
  }
  
  # calculate E and e
  temp_sum = matrix(0, nrow = p, ncol = p)
  temp_sum_scalar = 0
  for (i in 1:length(Xj))
  {
    Dinv = solve(Dj[[i]])
    XDinv = crossprod(Xj[[i]], Dinv)
    temp_sum = temp_sum + XDinv %*% Xj[[i]]
    temp_sum_scalar = temp_sum_scalar + XDinv %*% yj[[i]]
  }
  E = temp_sum + Binv
  e = temp_sum_scalar + Binv %*% a
  
  # calculate f
  temp_sum = 0
  for (i in 1:length(Xj))
  {
    temp_sum = temp_sum + crossprod(yj[[i]], solve(Dj[[i]])) %*% yj[[i]]
  }
  
  f = c(temp_sum - crossprod(e, solve(E)) %*% e + t(a) %*% solve(B) %*% a)
  
  top = -(N/2)*log(2*pi) + s*log(r) + lgamma(N/2 + s)
  bottom =  lgamma(s) + (c(determinant(x = B, logarithm = T)$modulus) + c(determinant(x = E, logarithm = T)$modulus) +
                           Reduce('+',lapply(X = Dj,FUN = function(X){c(determinant(x = X, logarithm = T)$modulus)})))/2
  inner = (N/2 + s)*log(f/2 + r)

  to_return = top - inner - bottom

  return(to_return)
}

evaluate_beta = function(beta, Cj)
{
  no_local = vector(length = length(Cj))
  for (i in 1:length(Cj))
  {
    temp = t(beta[[i]]) %*% solve(Cj[[i]]) %*% beta[[i]]
    no_local[i] = as.numeric(temp)
  }
  
  to_return = prod(no_local)
  return(to_return)
}

get_mNL_old = function(data, tau2 = 0.09/2)
{
  Xj = data$Xj
  yj = data$yj
  betaj = data$betaj
  nj = unlist(lapply(Xj, nrow))
  N = sum(nj)
  J = length(Xj)
  p = ncol(Xj[[1]])
  B = diag(p)
  Cj = get_Cj(Xj)
  a = rep(0, p)
  
  # some starting definitions 
  nu = N + 2*J
  
  # starting definitinos for G and g
  G = matrix(0, nrow = p, ncol = p)
  for (i in 1:J)
  {
    temp = crossprod(Xj[[i]])
    G = G + temp
  }
  G = G + solve(B)
  Ginv = solve(G)
  
  # redefinition for the stacked form (section 2.3)
  X = do.call(rbind, Xj)
  y = unlist(yj)
  XtX = lapply(Xj, crossprod) # list of X^TX
  XtX_whole = do.call(rbind, XtX) # stacked of XtX
  # Cj_inv = lapply(Cj, solve) # list of solve(C)
  Cj_inv = XtX
  
  # get M
  # start with getting the two diagonal poritons
  get_one = XtX[[1]]
  get_two = Cj_inv[[1]] / nj[1]
  for (i in 2:J)
  {
    temp_two = Cj_inv[[i]] / nj[i]
    get_one = adiag(get_one, XtX[[i]])
    get_two = adiag(get_two, temp_two)
  }
  get_two = (1/tau2) * get_two 
  
  M = get_one + get_two + (XtX_whole %*% Ginv %*% t(XtX_whole))  # TODO - fix this inverse 
  M = round(M, 5)
  Minv = round(solve(M), 5) # inverse of M, rounded for the precision 
  
  # get m
  Xty = crossprod((Xj[[1]]), yj[[1]])
  for (i in 2:J)
  {
    Xty = c(Xty, crossprod(Xj[[i]], yj[[i]]))
  }
  
  m = Xty + (XtX_whole %*% Ginv %*% (t(X) %*% y + solve(B) %*% a))
  
  # get form for t distribution
  mu = c(Minv %*% m)
  temp = (1/nu)*(sum(y^2) +
                   t(crossprod(X,y) + solve(B) %*% a) %*% 
                   Ginv %*% (crossprod(X,y) + solve(B) %*% a)
                 - (t(m) %*% Minv %*% m) + (t(a) %*% solve(B) %*% a))
  Epsilon = as.numeric(temp) * Minv
  # Epsilon = round(Epsilon, 5)
  
  
  
  top = 2^J * sterling_gamma(N/2 + J)
  C = prod(sqrt(unlist(lapply(Cj, det))))
  bottom = (pi^(N/2))*(p^J)*(tau2^(-J*(p/2 + 1)))*
    sqrt(det(G)*det(B)*det(M))*C
  const = (top/bottom) * (prod(nj))^(-(p/2 + 1))
  
  # sampler part
  values = 0
  max_rep = 1000
  for (i in 1:max_rep)
  {
    sample = c(LaplacesDemon::rmvt(n = 1, mu = mu, S = Epsilon, df = nu))
    sample = lapply(X = 1:J, FUN = function(X) {sample[((X-1)*p + 1):(X*p)]})
    values = values + evaluate_beta(sample, Cj)
  }
  
  int_part = values/max_rep
  BF = const*int_part
  to_return = as.numeric(log(BF))
  
  return(to_return)
}

get_mNL = function(data, tau2 = 0.09/2)
{
  Xj = data$Xj
  yj = data$yj
  betaj = data$betaj
  nj = unlist(lapply(Xj, nrow))
  N = sum(nj)
  J = length(Xj)
  p = ncol(Xj[[1]])
  B = diag(p)
  Cj = get_Cj(Xj)
  a = rep(0, p)
  
  # some starting definitions 
  nu = N + 2*J
  
  # starting definitinos for G and g
  G = matrix(0, nrow = p, ncol = p)
  g = matrix(0, nrow = ncol(Xj[[1]]), ncol = 1)
  yty = 0
  for (i in 1:J)
  {
    # getting G
    temp = crossprod(Xj[[i]])
    G = G + temp
    
    # getting g
    temp = crossprod(Xj[[i]], yj[[i]] - Xj[[i]] %*% betaj[[i]])
    g = g + g
  }
  G = G + solve(B)
  Ginv = solve(G)
  
  # redefinition for the stacked form (section 2.3)
  X = do.call(rbind, Xj)
  y = unlist(yj)

  # get M (equation 42)
  # start with getting the two diagonal poritons
  get_one = crossprod(Xj[[1]])
  get_two = solve(Cj[[1]]) / nj[1]
  XD = Xj[[1]]
  for (i in 2:J)
  {
    get_one = adiag(get_one, crossprod(Xj[[1]])) # getting diag(X_1^TX_1, ..., X_J^TX_J)
    get_two = adiag(get_two, solve(Cj[[i]]) / nj[i])
    XD = adiag(XD, Xj[[i]])
  }
  get_two = (1/tau2) * get_two 
  
  M = get_one + get_two - crossprod(XD, X) %*% Ginv %*% t(X) %*% XD # CHANGED TO MINUS
  M = round(M, 5)
  Minv = round(solve(M), 5) # inverse of M, rounded for the precision 
  
  # get m (equation 43)
  m = crossprod(XD, y) - crossprod(XD, X) %*% Ginv %*% (crossprod(X, y) + solve(B) %*% a) # CHANGED TO MINUS
  
  
  # get form for t distribution, (equation 45-46)
  mu = c(Minv %*% m)
  temp = crossprod(X, y) + solve(B) %*% a
  u = crossprod(y) - crossprod(temp, Ginv) %*% temp - crossprod(m, Minv) %*% m + crossprod(a, solve(B)) %*% a # CHNAGED TO MINUS
  Epsilon = as.numeric(u) * (1/nu) * Minv
  
  
  top = 2^J * sterling_gamma(nu/2) * tau2
  Cdet = prod(sqrt(unlist(lapply(Cj, det))))
  bottom = (pi^(N/2)) * (p^J) * u^(nu/2) * sqrt(det(G) * det(B) * det(M) * Cdet)
  const = (top/bottom) * (prod(nj))^(-(p/2 + 1))
  
  # sampler part
  values = 0
  max_rep = 1000
  for (i in 1:max_rep)
  {
    sample = c(LaplacesDemon::rmvt(n = 1, mu = mu, S = Epsilon, df = nu))
    sample = lapply(X = 1:J, FUN = function(X) {sample[((X - 1)*p + 1):(X*p)]})
    values = values + evaluate_beta(sample, Cj)
  }
  
  int_part = values/max_rep
  BF = const*int_part
  to_return = as.numeric(log(BF))
  
  return(to_return)
}

# bayes factors 
get_BF_NL = function(data, tau2 = 0.09/2)
{
  # get logs of each marginal
  null = get_m0(data, tau2 = 0.09/2)
  alt = get_mNL(data, tau2 = 0.09/2)
  
  # log BF is diffrence of non-local / local
  to_return = alt - null
  return(to_return)
}

get_BF_L = function(data, tau2 = 0.09/2)
{
  # get logs of each marginal
  null = get_m0(data, tau2 = 0.09/2)
  alt = get_mL(data, tau2 = 0.09/2)
  
  # log BF is diffrence of non-local / local
  to_return = null - alt
  return(to_return)
}

get_BF = function(data, tau2 = 0.09/2)
{
  # get logs of each marginal
  null = get_mL(data, tau2 = 0.09/2)
  alt = get_mNL(data, tau2 = 0.09/2)
  
  # log BF is diffrence of non-local / local
  to_return = null - alt
  return(to_return)
}

######
get_rho_comparison = function(prop_mean)
{
  rho = seq(0, 1, by = 0.1)
  BF = vector(length = length(rho))
  max_rep = 10
  for (i in 1:length(BF))
  {
    doParallel::registerDoParallel(cores = 80)
    tempBF = foreach::foreach(j = 1:max_rep, .combine = 'c', .multicombine = T) %dopar%{
      
      set.seed(j)
      temp = get_data(J = 10, rho = rho[i], prop_mean = prop_mean, sample_min = 5, sample_max = 5)
      
      get_BF_NL(temp)
    }
    
    BF[i] = mean(tempBF)
    print(rho[i])
  }
  return(list(rho = rho, BF =  BF))
}

tau2 = 0.3^2/2
set.seed(1)
p = 2
data = get_data(5, 0)
print(paste("value of local marginal: ", round(get_mL(data))))
print(paste("value of non-local marginal: ", round(get_mNL(data))))

print(paste("__________BF_________"))
print(get_BF(data))


# c = get_rho_comparison(0.3)
# plot(c$rho, c$BF, type = "l", main = paste("log(BF) verses correlation value", 0.5), ylab = "BF value", xlab = "rho")
