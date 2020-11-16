\name{picksamples}
\alias{picksamples}
\title{Pick subsamples to obtain initials}
\usage{
picksamples<-function(X,Y,n_ratio)
}
\description{
This function picks subsamples from the original data set (X,Y)
}
\arguments{
\item{X}{A data frame of explanatory variables. Intercept includes.}
\item{Y}{Response Variable.}
\item{n_ratio}{The ratio of data to be sampled.}
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
n_ratio=0.3
Sample<-picksamples(X,Y,n_ratio)
}
