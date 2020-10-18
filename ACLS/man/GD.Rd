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
p<-5
tau<-sqrt(n)/log(log(n))
mu<-matrix(0L,nrow=p,ncol=1)
#covarinace matrix
x_0<-matrix(1L,nrow=n,ncol=1)
rho<-0.5
Sig<-matrix(0L,nrow=p,ncol=p)
for (i in 1:p)
{
    for (j in 1:p){
    Sig[i,j]<-rho^abs(i-j)
}
  }
R<-mvrnorm(n, mu, Sig, tol = 1e-06, empirical = FALSE)
X_<-cbind(x_0,R)
X<-data.frame(X_)
colnames(X) <- c("Intercept", paste("X", 1:p, sep = ""))
beta_true<-c(0,3,4,1,2,0)
eps<-matrix(rnorm(n));
Y<-X_ %*% beta_true+eps
Y<-data.frame(Y)
beta_0<-c(0,0,0,0,0,0)
beta<-GD(beta_0,tau,X,Y)
}
