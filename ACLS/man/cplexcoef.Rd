\name{cplexcoef}
\alias{cplexcoef}
\title{Obtaining coefficients for CPLEX}
\usage{
cplexcoef<-function(X,Y,tau)
}
\description{
This function gets the coefficients for the MIQP problem to solve by Rcplex from data(X,Y)
}
\arguments{
\item{X}{A data frame of explanatory variables. Intercept includes.}
\item{Y}{Response Variable.}
\item{tau}{Adaptive Robustification Parameter. Could be a constant or a function of other parameter, for example, a function of sample size.}
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
coef<-cplexcoef(X,Y,tau)
}
