ipr<-function(X,y,print.it=FALSE,start=rep(1,dim(X)[2]),cutup=Inf,cutlow=cutup,epsrel=0.001,epsabs=0.1,maxiter=1000,det=FALSE){

  n<-dim(X)[1]
  p<-dim(X)[2]

  if(length(y)!=n) stop("Unable to use IPR. Length of y and number of rows of X differ.")
  if(length(start)!=p) stop("Unable to use IPR, Length of start and number of columns of X differ.")
  if(min(y)<=0) stop("Unable to use IPR. Values in y must be strictly positive.")
  if(min(X==0 | X==1)!=1) stop("Unable to use IPR. Matrix X must be filled with 0 and 1 only.")
  if(min(apply(X,1,max))==0) stop("Unable to use IPR. Some rows from matrix X have only 0's.")

  converge<-0
  I<-0

  mu_ipr<-start
  mu_ipr[mu_ipr<=0]<-min(mu_ipr[mu_ipr>0])/2

  while(converge==0){
    I<-I+1
    mutmp<-mu_ipr

    tmp1<-matrix(rep(y,each=p),n,p,byrow=T)
    tmp2<-matrix(rep(mu_ipr,n),n,p,byrow=T)
    Y<-X*tmp1*tmp2/matrix(rep(X%*%mu_ipr,p),n,p,byrow=F)

    if(I==1) Z<-matrix(1,n,p)

    if(I>=2){
      Blow<-matrix(rep(mu_ipr/cutlow,n),n,p,byrow=T)
      Bup<-matrix(rep(mu_ipr*cutup,n),n,p,byrow=T)
      Z<-1*(Y>=Blow)*(Y<=Bup)
    }

    mu_ipr<-apply(Z*Y,2,sum)/apply(Z*X,2,sum)

    epsrel0<-max(abs(mu_ipr-mutmp)/mutmp)
    epsabs0<-max(abs(mu_ipr-mutmp))

    converge<-1*(epsrel0<epsrel | epsabs0<epsabs | I==maxiter)

    if(print.it==TRUE) cat("ipr iteration",I,":",epsrel0,epsabs0,"\n")
  }

  tau <- apply(Y,2,sum)
  prop <- tau/sum(tau)

  if(det==TRUE){
    return(list(coef=mu_ipr,total=tau,proportions=prop,niter=I,epsrel=epsrel0,epsabs=epsabs0,detail=Y))
  } else {
    return(list(coef=mu_ipr,total=tau,proportions=prop,niter=I,epsrel=epsrel0,epsabs=epsabs0))
  }
}
