# A function to create a correlation matrix with a hierarchical structure

"make.hierarchical" <-
function (gload=NULL,fload=NULL,n=0,raw=FALSE) {
 require(MASS)
 if(is.null(gload)) gload=matrix(c(.9,.8,.7),nrow=3)
 if(is.null(fload)) {fload <-matrix(c(               
            .8,0,0,
            .7,0,.0,
            .6,0,.0,
             0,.7,.0,
             0,.6,.0,
             0,.5,0,
              0,0,.6,
              0,0,.5,
              0,0,.4),   ncol=3,byrow=TRUE)}
  fcor <- gload %*% t(gload)           #the factor correlation matrix
  diag(fcor) <-1                       #put ones on the diagonal
  model <-  fload %*% fcor %*% t(fload) #the model correlation matrix for oblique factors
  diag(model)<- 1                       # put ones along the diagonal 
  if(n>0) {
  	model <- mvrnorm(n = n, Sigma=model, tol = 1e-6, empirical = FALSE)
  	if (!raw ) { model <- cor(model) } }
  make.hierarchical <- model }
