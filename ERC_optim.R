#' ---
#' title: "Equally-Weighted Risk Contributions by optimization "
#' author: "Richard Warnung"
#' date: "August 17, 2022"
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

#' # Equally-Weighted Risk Contributions (ERC) portfolio
#' Optimization with `constrOptim.nl` from package `alabama`.
#' We first solve the problem
#' $$\sqrt{y^T \Sigma y} \rightarrow Min$$ 
#' subject to
#' \begin{align*}
#'    y_i > 0 \\
#'    \sum \ln y_i \ge c
#'\end{align*} 
#' Then we rescale $x_i = y_i/\sum_{i=1}^n y_i$ to get the ERC portfolio.
#' 
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

#' ## Conclusion
#' 
#' As long as $c$ is below the upper bound $-n \log(n)$, solutions $y$ are different but lead to the ERC after scaling.
#' 
#' # Direct derivation of the ERC portfolio without scaling
#' This time we solve the following optimization problem: 
#' $$\sqrt{x^T \Sigma x} \rightarrow Min$$ 
#' subject to
#' \begin{align*}
#'    x_i &> 0 \\
#'    \sum_{i=1}^n x_i &= 1 \\
#'    \sum_{i=1}^n \ln x_i &\ge c^*.
#' \end{align*}
#' 
#' We define a new constraint function and fix the constant to $c^{\ast} = c - n \ln( \sum_{i=1}^n y_i^{\ast})$. Where $c$ was used to derive the optimal solution $y^{\ast}$.
#' According to Teiletche and Roncalli this will result in the ECR directly without the need to rescale.
#+ 
c_star = c - length(y3)*log(sum(y3))
constr_fn_star = function(par){
  h = rep(NA, length(par)+1)
  h[1:length(par)] = par # par[i] > 0
  h[length(par)+1] =  sum(log(par))-c_star  # sum(log(par)) > c*
  return( h )
}

constr_eq_fn = function(par){
  h = rep(NA,1)
  h[1] = sum(par)-1 
}

ans = constrOptim.nl(par = c(0.18, 0.22, 0.29, 0.31), fn = objective_fn, hin = constr_fn_star, heq = constr_eq_fn, control.outer = list(trace = FALSE))  

x = ans$par
x
sum(x)

vol_ERC = sqrt( t(x)%*%sigma.mat%*%x)[1,1]
vol_ERC
contribs = x *(sigma.mat%*%x/vol_ERC)
contribs


#' We see that the solution $x$ gives the ERC (apart from numerical inaccuracies).
#' 
#' ## Changing the constant $c^*$
#' 
#' Next we investigate what happens if we change the constant c_star. We go back to the initial constant $c^* = c$.
#+

c_star = c 
ans = constrOptim.nl(par = c(0.18, 0.22, 0.29, 0.31), fn = objective_fn, hin = constr_fn_star, heq = constr_eq_fn, control.outer = list(trace = FALSE))  

x = ans$par
x
sum(x)

vol_pf = sqrt( t(x)%*%sigma.mat%*%x)[1,1]
vol_pf
contribs = x *(sigma.mat%*%x/vol_pf)
contribs

#' We see that the portfolio constraint is still fulfilled, but the ERC-property is not fulfilled. The volatility is lower than in the case of the ERC.
#' 
#' ## Removing the portfolio constraint
#' 
#' Next we investigate the result if we remove the portfolio constraint $\sum x_i = 1$:

c_star = c - length(y3)*log(sum(y3))
ans = constrOptim.nl(par = c(0.18, 0.22, 0.29, 0.31), fn = objective_fn, hin = constr_fn_star, control.outer = list(trace = FALSE))  

x = ans$par
x
sum(x)

vol_pf = sqrt( t(x)%*%sigma.mat%*%x)[1,1]
vol_pf
contribs = x *(sigma.mat%*%x/vol_pf)
contribs

#' The ERC and the portfolio property are fulfilled.
#' 
#' 

#' ## Changing the objective function
#' 
#' Due to the dominance of the condition $\sum_{i=1}^n \ln x_i \ge c^{\ast}$ we investigate whether we can replace the objective function by another convex function. We choose
#' the L2-norm $\sum_{i=1}^n x_i^2$.
#' 
#+
#+ 

objective_fn = function(par){
  return(  sum(par^2) ) ## norm as objective function
}

ans = constrOptim.nl(par = c(0.18, 0.22, 0.29, 0.31), fn = objective_fn, hin = constr_fn_star, control.outer = list(trace = FALSE))  

x = ans$par
x
sum(x)

vol_pf = sqrt( t(x)%*%sigma.mat%*%x)[1,1]
vol_pf
contribs = x *(sigma.mat%*%x/vol_pf)
contribs

#' Changing the objective function destroys the ERC property.
#'
#' ## Derivation of the mathematics
#' 
#' Let us see why $c^{\ast} = c - n \ln( \sum_{i=1}^n y_i^{\ast})$ gives the correct solution. This part of the post roughly reflects personal communications with in the of
#' the authors of the mentioned article.
#' 
#' We consider $x^* := y^* \phi$, where $\phi = 1/\sum_{i=1}^n y^*_i$. Thus $x^*$ is the re-scaled ERC-portfolio satisfying $\sum_{i=1}^n x^*_i = 1$.
#' 
#' By definition of $x^*$, it holds that
#' $$  \sqrt{x^{*T} \Sigma x^*} =  \phi \sqrt{y^{*T} \Sigma y^*}.$$ 
#' Thus, the resulting volatility is scaled by $\phi$.
#' As 
#' $$y_i^{\ast} > 0$$ for $i=1,\ldots,n$, also $\phi>0$ and thus $$x^{*}_i > 0$$ for $i=1,\ldots,n$.
#' Finally, looking at the logarithmic constraint we have from the optimality of $y^*$ that
#' \begin{align*}
#'  \sum_{i=1}^n \ln y_i^* &\ge c \\
#'  \sum_{i=1}^n \ln \left(x_i^*\phi \right) &\ge c \\
#'  \sum_{i=1}^n \ln(x_i^*) - \sum_{i=1}^n \ln( \phi) &\ge c \\
#'  \sum_{i=1}^n \ln(x_i^*) &\ge c + n \ln( \phi)
#'\end{align*}
#' and as $\phi = 1/\sum_{i=1}^n(y^*_i)$ we arrive at
#' $$
#' \sum_{i=1}^n \ln(x_i^*) \ge c - n \ln(\sum_{i=1}^n y^*_i).
#' $$


#+ include=FALSE
options(warn = oldw)