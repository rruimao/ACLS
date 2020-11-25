# ACLS
Adaptive Capped Least Squares
## Description
This package includes two methods applied to minimize the adaptive capped least squares loss: randomized gradient descent method and gradient descent method with initials obtained from CPLEX.

Suppose we observe data vectors  (x<sub>i</sub>,y<sub>i</sub>) that follow a linear model y<sub>i</sub>=x<sub>i</sub><sup>T</sup>&beta;<sup>*</sup>+&epsilon;<sub>i</sub>, i=1,...n, where y<sub>i</sub> is a univariate response,  x<sub>i</sub> is a d-dimensional predictor, &beta;<sup>*</sup> denotes the vector of regression coefficients, and &epsilon;<sub>i</sub> is a random error. We propose the adpative capped least squares loss, &ell;(x)=x<sup>2</sup>/2 if |x| &leq; &tau;; &tau;<sup>2</sup>/2, if |x| &gt; &tau;, where &tau;=&tau;(n) &gt; 0 is referred to as the adaptive capped least squares parameter. The proposed methods are applied to find &beta; that minimizes L(&beta;)= n<sup>-1</sup> &sum; l(y<sub>i</sub>-x<sub>i</sub><sup>T</sup> &beta; ).

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
We present two examples: random generated data with y-outliers and random generated data with x-outliers and y-outliers. 



### First example: random generated data with y-outliers
we generate contaminated random errors &epsilon;<sub>i</sub> from a mixture of normal distribution 0.9N(0,1)+0.1N(10,1) and x<sub>i</sub>'s are independently and identically distributed (i.i.d.) from N(0,&Sigma;) where &Sigma;=0.5<sup>|j-k|</sup>. We set &beta;<sup>*</sup> =(0,3,4,1,2,0)<sup>T</sup> to generate y<sub>i</sub>. We provide one example of this type, "ex_1.Rdata", and it can be downloaded from example file.

```R
load("ex_1.Rdata")
Y<-ex_1$Y
X<-ex_1[,!(names(ex_1) %in% "Y")]
```

We first apply our proposed methods, randomized gradient descent (denoted by ACLS) and gradient descent with intials obtained from CPLEX (denoted by ACLS-h), to the data. 


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
We then compare mean square errors (MSEs), defined as || &beta;&circ;-&beta;<sup>*</sup>||<sup>2</sup>, of these two methods with MSEs of ordinary least squares method (OLS), Huber method with adaptive resistant parameter (denoted by AHR) and least trimmed squares method (LTS). Results obtained from the CPLEX on the whole dataset (denoted by ACLS-C) are treated as benchmark. 
``` R
install.packages("robustbase")
library("robustbase")
install.packages("robustreg")
library("robustreg")
Z.OLS<-lm(Y~X1+X2+X3+X4+X5,data=ex_1)
beta.OLS<-Z.OLS$coefficients
Z.Huber<-robustRegH(Y~X1+X2+X3+X4+X5,data=ex_1,tune=tau)
beta.Huber<-Z.Huber$coefficients
Z.LTS<-ltsReg(X[,-1],Y,intercept=TRUE,adjust=TRUE)
beta.LTS<-Z.LTS$coefficients
Coef<-cplexcoef(X,Y,tau)
ans<-Rcplex(Coef$f,Coef$A,Coef$b,Coef$Q,Coef$lb,Coef$ub,vtype=Coef$vtype)
beta.cplex<-ans$xopt[(3*n+1):(3*n+d+1)]
```

We summarize the MSEs of all methods in the following table.

|    |OLS | AHR |  LTS | ACLS | ACLS-h | ACLS-C |
| :---         |     :---:      |          ---: |          ---: |          ---: |          ---: |          ---: |
| MSE   | 4.3016    | 4.3016 |0.5026|0.0626|0.0625|0.0636|

### Second example: random generated data with $x$-outliers and $y$-outliers
we generate contaminated random errors &epsilon;<sub>i</sub> from a mixture of normal distribution 0.9N(0,1)+0.1N(10,1) and x<sub>i</sub>'s are independently and identically distributed (i.i.d.) from N(0,&Sigma;) where &Sigma;=0.5<sup>|j-k|</sup>. We then add a random perturbation vector  z<sub>i</sub> ~ N(10 &times; 1<sub>d-1</sub>,I<sub>d-1</sub>) to each covariate x<sub>i</sub> in the contaminated samples. We also use &beta;<sup>*</sup> =(0,3,4,1,2,0)<sup>T</sup> and use uncontaminated x<sub>i</sub> to generate y<sub>i</sub>. We provide one example of this type, "ex_2.Rdata", and it can be downloaded from example file.
	
``` R
load("ex_2.Rdata")
Y<-ex_2$Y
X<-ex_2[,!(names(ex_1) %in% "Y")]
```	


We use the same code in the first example with the new data to get estimators for all methods. We also collect the MSEs in the following table. 

|    |OLS | AHR |  LTS | ACLS | ACLS-h | ACLS-C |
| :---         |     :---:      |          ---: |          ---: |          ---: |          ---: |          ---: |
| MSE   | 18.6446   | 18.6446    |0.7033|0.3558|0.3555|0.3569|

## Reference
Huber, P. J. (1973). Robust regression: asymptotics, conjectures and Monte Carlo. *The Annals of Statistics* **1** 799-821.

Rousseeuw, P. J. and Driessen, K. V. (1999). A fast algorithm for the minimum covariance determinant estimator. *Technometrics* **41** 212-223.

Rousseeuw, P. J. (1984). Least median of squares regression. *Journal of the American Statistical Association* **79** 871-880.

Sun, Q., Zhou, W.-X. and Fan, J. (2020). Adaptive Huber regression. *Journal of the American Statistical Association* **115** 254-265.

