#Written by Andreas Moltner with some revisions by William Revelle
#March, 2010
"glb.algebraic"<- 
function(Cov,LoBounds=NULL, UpBounds=NULL)
{ if(!requireNamespace('Rcsdp')) {stop("Rcsdp must be installed to find the glb.algebraic")
     }
cl<-match.call()
  # check input
  p<-dim(Cov)[2]
 if (dim(Cov)[1] != p) Cov <- cov(Cov) #find the covariances
  if (any(t(Cov)!=Cov))
    stop("'Cov' is not symmetric")
   if(is.null(LoBounds)) LoBounds <-rep(0,ncol(Cov))
   if(is.null(UpBounds)) UpBounds <- diag(Cov) 
  
  if (any(LoBounds>UpBounds))
  { stop("'LoBounds'<='UpBounds' violated")
  }
#  if (min(eigen(Cov,symmetric=TRUE,only.values=TRUE)$values)<0)
#    stop("'Cov' is not positive semidefinite")
  if (length(LoBounds) != p) 
    stop("length(LoBounds) != dim(Cov)")
  if (length(UpBounds) != p) 
    stop("length(UpBounds)!=dim(Cov)")
  Var<-diag(Cov)
  # objective function opt --> min
  opt=rep(1,p)
  # set up csdp input
  C<-list(diag(Var)-Cov,
          -UpBounds,LoBounds)
  A<-vector("list",p)
  for (i in 1:p)
  { b<-rep(0,p)
    b[i]<-1
    A[[i]]<-list(diag(b),-b,b)
  }
  K<-list(type=c("s","l","l"),size=rep(p,3))
  # call csdp
  result<- Rcsdp::csdp(C,A,opt,K)
  if (result$status>=4||result$status==2)
  { warning("Failure of csdp, status of solution=",result$status)
    lb<-list(glb=NA,solution=NA,status=result$status,Call=cl)
  } else
  { if (result$status!=0)
    { warning("status of solution=",result$status)
    }
    # greatest lower bound to reliability
    item.diag <- result$y
    names(item.diag) <- colnames(Cov)
    lb<-list(glb=(sum(Cov)-sum(Var)+sum(result$y))/sum(Cov),
            solution = item.diag,
            status=result$status,
            Call=cl)
  }
  return(lb)
}  



