picksamples<-function(X,Y,n_ratio){
  n<-nrow(X)
  p<-ncol(X)-1
  X.matrix<-data.matrix(X, rownames.force = NA)
  Y.matrix<-data.matrix(Y, rownames.force = NA)
  n_sp<-round(n*n_ratio)
  X_sp<-matrix(0L,nrow=round(n_sp),ncol=p+1)
  Y_sp<-matrix(0L,nrow=round(n_sp),ncol=1)
  index<-matrix(sample(1:n,round(n_sp)))
  for (i in 1:dim(index)[1])
  {
    j<-index[i]
    X_sp[i,]<-X.matrix[j,]
    Y_sp[i,]<-Y.matrix[j,]
  }
  Sample<-list();
  Sample[[1]]<-X_sp;
  Sample[[2]]<-Y_sp;
  names(Sample)<-c("X_sp","Y_sp");
  return(Sample);
}
