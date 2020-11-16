#generate initials
#to compare the loss and choose the min loss as initial
RGD<-function(X,Y,tau,iter,eta_0=1e-3,alpha=2){
  library("uniformly")
  n<-nrow(X)
  p<-ncol(X)-1
  X.matrix<-data.matrix(X, rownames.force = NA)
  Y.matrix<-data.matrix(Y, rownames.force = NA)
  beta_gdm<-matrix(0L,nrow=iter,ncol=p+1)
  L<-matrix(0L,nrow=iter,ncol=1)
  for (s in 1:iter){
  beta_0=t(runif_in_pball(1,p+1,2,tau))
  beta_gdm[s,]<-GD(beta_0,tau,X.matrix,Y.matrix,eta_0,alpha)
    r<-Y.matrix-X.matrix%*%beta_gdm[s,]
    L[s]<-1/2*sum((r[which(abs(r)<tau)])^2)+1/2*tau^2*length(r[which(abs(r)>tau)])
  }
  beta<-matrix(beta_gdm[which(L==min(L)),])
  beta1<-beta[,1]
  beta_fin<-GD(beta_0,tau,X.matrix,Y.matrix,eta_0,alpha)
  return(beta_fin)
}







