#' ---
#' title: "ERC by optimization "
#' author: "Richard Warnung"
#' date: "June 6, 2022"
#'---

#' # Introduction
#' The approach is based on [On the Properties of Equally-Weighted Risk Contributions Portfolios](https://www.researchgate.net/publication/45397778_On_the_Properties_of_Equally-Weighted_Risk_Contributions_Portfolios)
#' by Teiletche and Roncalli. \
#' Required packages: alabama 
#'
#' # Generating data
#' 4 assets with decreasing volatility and equal correlation:
#+

vols = matrix( c(0.3, 0.25, 0.2, 0.15), nrow = 4)

rho = 0.7
cor.mat = matrix(c(rho), nrow=4, ncol=4)
diag(cor.mat) = 1

sigma.mat =  cor.mat * matrix(outer(vols, vols), nrow=4, byrow=TRUE) 

sigma.mat
## checking the covariance matrix
cov2cor(sigma.mat)

#' # Equally weighted portfolio
#+
x_equal = matrix(c(1/4, 1/4, 1/4, 1/4), nrow=4)

#' volatility and contributions of equally weighted portfolio:
#+ 
vol_equal = sqrt( t(x_equal)%*%sigma.mat%*%x_equal )[1,1]
vol_equal
contribs_equal_w = x_equal *(sigma.mat%*%x_equal/vol_equal)
contribs_equal_w

#' # Equal Risk Contributions (ERC) portfolio
#' Optimization with `constrOptim.nl` from package `alabama`.

#+ 
library(alabama)

objective_fn = function(par){
  return( (t(par)%*%sigma.mat%*%par)[1,1]) ## no square-root, same arg min
}
c_max = -nrow(sigma.mat)*log(nrow(sigma.mat))

constr_fn = function(par){
  h = rep(NA, length(par)+1)
  h[1:length(par)] = par # par[i] > 0
  h[length(par)+1] =  sum(log(par))-c  # sum(log(par)) > c
  return( h )
}

oldw = getOption("warn") ## silencing the warnings during optimization, keeping the settings for later.
options(warn = -1)

#' ## Choosing c = c_max-0.4
#+
c = c_max-0.4

ans = constrOptim.nl(par = c(1/5, 1/5, 1/5, 1-3/5), fn = objective_fn, hin = constr_fn, control.outer = list(trace = FALSE))  
y = ans$par
y
sum(y)
x = y/sum(y)
x
vol_ERC = sqrt( t(x)%*%sigma.mat%*%x)[1,1]
vol_ERC

contribs = x *(sigma.mat%*%x/vol_ERC)
contribs

#' ## Choosing c = c_max-0.5
#+
c = c_max-0.5

ans = constrOptim.nl(par = c(1/5, 1/5, 1/5, 1-3/5), fn = objective_fn, hin = constr_fn,  control.outer = list(trace = FALSE))  
y2 = ans$par
y2
sum(y2)
x = y2/sum(y2)
x
vol_ERC = sqrt( t(x)%*%sigma.mat%*%x)[1,1]
vol_ERC
contribs = x *(sigma.mat%*%x/vol_ERC)
contribs

#' ## Choosing c = c_max-1
#+
c = c_max-1

ans = constrOptim.nl(par = c(1/5, 1/5, 1/5, 1-3/5), fn = objective_fn, hin = constr_fn,  control.outer = list(trace = FALSE))  
y3 = ans$par
y3
sum(y3)
x = y3/sum(y3)
x
vol_ERC = sqrt( t(x)%*%sigma.mat%*%x)[1,1]
vol_ERC
contribs = x *(sigma.mat%*%x/vol_ERC)
contribs

options(warn = oldw)
#' ## Conclusion
#' As long as $c$ is below the upper bound $-n \log(n)$, solutions $y$ are different but lead to the ERC after scaling.
