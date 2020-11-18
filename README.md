# ACLS
Adaptive Capped Least Squares
## Description
This package includes two methods applied to minimize the adaptive capped least squares loss: randomized gradient descent method and gradient descent method with initials obtained from CPLEX.

Suppose we observe data vectors  (x<sub>i</sub>,y<sub>i</sub>) that follow a linear model y<sub>i</sub>=x<sub>i</sub><sup>T</sup>&beta;<sup>*</sup>+&epsilon;<sub>i</sub>, i=1,...n, where y<sub>i</sub> is a univariate response,  x<sub>i</sub> is a d-dimensional predictor, &beta;<sup>*</sup> denotes the vector of regression coefficients, and &epsilon;<sub>i</sub> is a random error. We propose the adpative cappled least squares loss, l(x)=x<sup>2</sup>/2 if |x| &leq; &tau; &tau;<sup>2</sup>/2, if |x| &gt; &tau; where &tau;=&tau;(n) &gt; 0 is referred to as the adaptive capped least squares parameter. The proposed methods are applied to find &beta; that minimizes L(&beta;)= n<sup>-1</sup> &sum; l(y<sub>i</sub>-x<sub>i</sub><sup>T</sup> &beta; ).

## Installation
Install **ACLS** from GitHub:
``` R
install.packages("devtools")
library(devtools)
devtools::install_github("rruimao/ACLS/ACLS")
library(ACLS)
``` 

Install **Rcplex** if the user wants to use CPLEX to obtain the intial or the solution to the problem:

See the **INSTALL** file from https://cran.r-project.org/web/packages/Rcplex/index.html.


## Functions
There are four functions in this package:

- **GD**: Gradient descent method with user-defined initial.
- **RGD**: Randomized gradient descent method.
- **picksamples**: Pick subsamples to obtain the initial.
- **cplexcoef**: Get the coefficients for the problem solved by CPLEX.

## Examples
We present two examples: random generated data with y-outliers and random generated data with x-outliers and x-outliers. 



### First example: random generated data with y-outliers
we generate contaminated random errors &epsilon;<sub>i</sub> from a mixture of normal distribution 0.9 N(0,1)+0.1N(10,1) and x<sub>i</sub>'s are independently and identically distributed (i.i.d.) from N(0,&Sigma;) where &Sigma;=0.5<sup>|j-k|</sup>. We set &beta;<sup>*</sup> =(0,3,4,1,2,0)<sup>T</sup> to generate y<sub>i</sub>. This random genrated data can be downloaded from example file.

```R
load
```

We first apply our proposed methods, randomized gradient descent (denoted by ACLS) and gradient descent with intials obtained from CPLEX (denoted by ACLS-h), to the data. We then compare mean square errors (MSEs) of these two methods with MSEs of ordinary least squares method (OLS), Huber method with adaptive resistant parameter (denoted by AHR) and least trimmed squares method (LTS). Results obtained from the CPLEX on the whole dataset (denoted by ACLS-C) are treated as benchmark. 


### Randomized gradient descent
We randomly generate 10 initials &beta;<sup>*</sup> ~ Unif(B<sub>2</sub>(&tau;)), where Unif(B<sub>2</sub>(&tau;)) is a uniform distribution on the l<sub>2</sub>-ball B<sub>2</sub>(&tau;)={x: ||x||<sub>2</sub> &leq; &tau; }. This method finds the initial that provides the smallest adaptive capped least squares loss.

``` R
tau<-sqrt(n)/log(log(n))
iter<-10
beta.rgd<-RGD(X,Y,tau,iter)
```

### Gradient descent method with initials obtained from CPLEX 
We first sample 30% data and use CPLEX to find the optimal solution as the initial, then we apply gradient descent method with this initial.
``` R
library("Rcplex")
n_ratio<-0.3
n_sp<-round(n*n_ratio)
Sample<-picksamples(X,Y,n_ratio)
Coef_sp<-cplexcoef(Sample$X_sp,Sample$Y_sp,tau)
ans_sp<-Rcplex(Coef_sp$f,Coef_sp$A,Coef_sp$b,Coef_sp$Q,Coef_sp$lb,Coef_sp$ub,vtype=Coef_sp$vtype)
beta_0<-ans_sp$xopt[(3*n_sp+1):(3*n_sp+d+1)]
beta.gdc<-GD(beta_0,tau,X,Y)
```

### Other methods
``` R
install.packages("robustbase")
library("robustbase")
install.packages("robustreg")
library("robustreg")
beta.OLS<-solve(t(X)%*%X)%*%t(X)%*%Y
Z<-data.frame(X,Y)
colnames(Z)<-c("X1","X2","X3","X4","X5","X6","Y")
Z.Huber<-robustRegH(Y~X2+X3+X4+X5+X6,data=Z,tune=tau)
beta.Huber<-Z.Huber$coefficients
Z.LTS<-ltsReg(X[,-1],Y,intercept=TRUE,adjust=TRUE)
beta.LTS<-Z.LTS$coefficients
Coef<-cplexcoef(X,Y,tau)
ans<-Rcplex(Coef$f,Coef$A,Coef$b,Coef$Q,Coef$lb,Coef$ub,vtype=Coef$vtype)
beta.cplex<-ans$xopt[(3*n+1):(3*n+d+1)]
```

We summerize the MSEs of all methods in the following table.

|    |OLS | AHR |  LTS | ACLS | ACLS-h | ACLS-C |
| :---         |     :---:      |          ---: |          ---: |          ---: |          ---: |          ---: |
| MSE   | 4.3016    | 4.3016 |0.5026|0.0626|0.0625|0.0636|

### Second example: random generated data with $x$-outliers and $y$-outliers
we generate contaminated random errors &epsilon;<sub>i</sub> from a mixture of normal distribution 0.9 N(0,1)+0.1N(10,1)$ and x<sub>i</sub>'s are independently and identically distributed (i.i.d.) from N(0,&Sigma;) where &Sigma;=0.5<sup>|j-k|</sup>. We then add a random perturbation vector  z<sub>i</sub> ~ N(10 &times; 1<sub>d-1</sub>,I<sub>d-1</sub>) to each covariate x<sub>i</sub> in the contaminated samples. We also use &beta;<sup>*</sup> =(0,3,4,1,2,0)<sup>T</sup> and use uncontaminated x<sub>i</sub> to generate y<sub>i</sub>.
	
``` R
load
```	


We use the same code in the first example replacing X, Y with X_new, Y_2 respectively to get estimators for all methods. We also collect the MSEs in the following table.

|    |OLS | AHR |  LTS | ACLS | ACLS-h | ACLS-C |
| :---         |     :---:      |          ---: |          ---: |          ---: |          ---: |          ---: |
| MSE   | 18.6446   | 18.6446    |0.7033|0.3558|0.3555|0.3569|
