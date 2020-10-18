#gradient descent method with user-defined intial
GD<-function(beta_0,tau,X,Y,eta_0=1e-3,alpha=2){
  n<-nrow(X)
  p<-ncol(X)-1
  X.matrix<-data.matrix(X, rownames.force = NA)
  Y.matrix<-data.matrix(Y, rownames.force = NA)

  e_0<-Y.matrix-X.matrix%*%beta_0
  L_0<-(1/2*sum((e_0[which(abs(e_0)<tau)])^2)+1/2*tau^2*length(e_0[which(abs(e_0)>tau)]))/n
  w<-matrix(1L,nrow=n,ncol=1)
  w[which(e_0<(-tau))]<-0
  w[which(e_0>tau)]<-0
  ew<-e_0*w
  diff<--t(X.matrix)%*%ew/n

  beta<-beta_0-eta_0*diff
  e_1<-Y.matrix-X.matrix%*%beta
  L_1<-(1/2*sum((e_1[which(abs(e_1)<tau)])^2)+1/2*tau^2*length(e_1[which(abs(e_1)>tau)]))/n

  eta<-alpha*eta_0

  beta<-beta_0-eta*diff
  e_2<-Y.matrix-X.matrix%*%beta
  L_2<-(1/2*sum((e_2[which(abs(e_2)<tau)])^2)+1/2*tau^2*length(e_2[which(abs(e_2)>tau)]))/n

    while (L_2<L_1){
      L_1<-L_2
      eta<-alpha*eta
      beta<-beta_0-eta*diff
      e_2<-Y.matrix-X.matrix%*%beta
      L_2<-(1/2*sum((e_2[which(abs(e_2)<tau)])^2)+1/2*tau^2*length(e_2[which(abs(e_2)>tau)]))/n
      }
    eta=eta/alpha;
    beta<-beta_0-eta*diff

  if (L_0<L_1){
    beta<-beta_0
    }

    j<-1

    while ((t(beta-beta_0)%*%(beta-beta_0)>1e-8 && j<10000000)){
      beta_0=beta
      e_0<-Y.matrix-X.matrix%*%beta_0
      L_0<-(1/2*sum((e_0[which(abs(e_0)<tau)])^2)+1/2*tau^2*length(e_0[which(abs(e_0)>tau)]))/n
      w<-matrix(1L,nrow=n,ncol=1)
      w[which(e_0<(-tau))]<-0
      w[which(e_0>tau)]<-0
      ew<-e_0*w
      diff<--t(X.matrix)%*%ew/n
      beta<-beta_0-eta_0*diff
      e_1<-Y.matrix-X.matrix%*%beta
      L_1<-(1/2*sum((e_1[which(abs(e_1)<tau)])^2)+1/2*tau^2*length(e_1[which(abs(e_1)>tau)]))/n
      eta<-alpha*eta_0
      beta<-beta_0-eta*diff
      e_2<-Y.matrix-X.matrix%*%beta
      L_2<-(1/2*sum((e_2[which(abs(e_2)<tau)])^2)+1/2*tau^2*length(e_2[which(abs(e_2)>tau)]))/n

      while (L_2<L_1){
          L_1<-L_2
          eta<-alpha*eta
          beta<-beta_0-eta*diff
          e_2<-Y.matrix-X.matrix%*%beta
          L_2<-(1/2*sum((e_2[which(abs(e_2)<tau)])^2)+1/2*tau^2*length(e_2[which(abs(e_2)>tau)]))/n
          }
        eta=eta/alpha;
        beta<-beta_0-eta*diff
        if (L_0<L_1){
          beta<-beta_0
        }
        j<-j+1
    }





  return(beta)
}
