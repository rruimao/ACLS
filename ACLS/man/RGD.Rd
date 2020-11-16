\name{RGD}
\alias{RGD}
\title{Randomized Gradient Descent Method}
\usage{
RGD(X,Y,tau,iter,eta_0=1e-3,alpha=2)
}
\description{
Gradient descent method with random initial (generated within d-ball) that provides smallest adaptive resistant loss.
}
\arguments{
\item{X}{A data frame of explanatory variables. Intercept includes.}
\item{Y}{Response Variable}
\item{tau}{Adaptive Robustification Parameter. Could be a constant or a function of other parameter, for example, a function of sample size}
\item{iter}{The number of iterations to choose the initial}
\item{eta_0}{The initial step size}
\item{alpha}{Inflation factor}
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
iter<-10
beta<-RGD(X,Y,tau,10)
}
