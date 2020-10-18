#get the coefficients for the MIQP problem to solve by Rcplex from data(X,Y)
cplexcoef<-function(X,Y,tau){
  n<-nrow(X)
  p<-ncol(X)-1
  X.matrix<-data.matrix(X, rownames.force = NA)
  Y.matrix<-data.matrix(Y, rownames.force = NA)
  Q=matrix(0L,nrow=3*n+p+1,ncol=3*n+p+1)
  for (i in 1:n)
  {
    Q[i,i]<-2
  }

  f=matrix(0L,nrow=1,ncol=3*n+p+1)
  f[,(2*n+1):(3*n)]=tau^2
  M=10000

  I=diag(n)
  Z=matrix(0L,nrow=n,ncol=n)
  Z1=matrix(0L,nrow=n,ncol=p+1)
  A1=cbind(-I,-I,Z,-X.matrix)
  A2=cbind(-I,-I,Z,X.matrix)
  A3=cbind(Z,I,-M*I,Z1)
  A=rbind(A1,A2,A3)
  Z2=matrix(0L,nrow=n,ncol=1)
  b=rbind(-Y.matrix,Y.matrix,Z2)

  O=matrix(1L,nrow=n+p+1,ncol=1)
  O1=matrix(1L,n,1)
  O2=matrix(1L,2*n+p+1,1)
  lb=rbind(Z2,Z2,-Inf*O)
  ub=rbind(tau*O1,Inf*O2)

  vtype=""
  for (i in 1:(2*n)){
    vtype=rbind(vtype,"C")
  }
  for (i in (2*n+1):(3*n)){
    vtype=rbind(vtype,"B")
  }
  for (i in (3*n+1):(3*n+p+1)){
    vtype=rbind(vtype,"C")
  }
  vtype=vtype[-1]
  vtype=t(vtype)

  Coef<-list();
  Coef[[1]]<-Q;
  Coef[[2]]<-f;
  Coef[[3]]<-A;
  Coef[[4]]<-b;
  Coef[[5]]<-lb;
  Coef[[6]]<-ub;
  Coef[[7]]<-vtype
  names(Coef)<-c("Q","f","A","b","lb","ub","vtype");
  return(Coef);
}
