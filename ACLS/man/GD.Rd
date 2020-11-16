\name{GD}
\alias{GD}
\title{Steepest Gradient Descent Approximation Method}
\usage{
GD<-function(beta_0,tau,X,Y)
}
\description{
The steepest gradient descent approximation method find the approximation of the exact step size to minimize the loss function, and use gradient descent method to obtain the estimate
}
\arguments{
\item{beta_0}{The initial of the coefficient.}
\item{X}{A data frame of explanatory variables. Intercept includes.}
\item{Y}{Response Variable.}
\item{tau}{Adaptive Robustification Parameter. Could be a constant or a function of other parameter, for example, a function of sample size.}
\item{eta_0}{The initial step size.}
\item{alpha}{Inflation factor.}
}
\examples{
n<-50
d<-5
a=matrix(data = rnorm(1, 0, 1), nrow = n, ncol = d)
x_0<-matrix(1L,nrow=n,ncol=1)
X=cbind(x_0,a)
beta_true<-c(0,3,4,1,2,0)
eps<-matrix(rnorm(n));
#Genarate response Y using true coefficient beta_true
Y<-X %*% beta_true+eps
tau<-sqrt(n)/log(log(n))
beta_0<-c(0,rnorm(5))
beta<-GD(beta_0,tau,X,Y)
}
