

#Developed 6/15/21
#Fixed 6/20/21 to avoid the problem of scales of one item or too many to find splits
#Modified 7/10/21 to allow it to function with just a set of items, no keys specified
#Modified 7/11/21 to include the unidimensional estimates from unidim  and to use the R object for speed

#Further modified 9/30/23 to force 1 factor solution length(keys < 4)
#added the covar option 10/23/23
#added fm option  2/29/24
reliability <- function(keys=NULL,items,nfactors=2,n.obs=NA,split=TRUE,raw=TRUE,plot=FALSE,hist=FALSE,
   n.sample=10000,brute=FALSE,check.keys=TRUE,covar=FALSE,fm="minres")  {
   cl <- match.call()
  result <- list()
  splits <- list()
  best <- worst <- list()
   res.name <- list()
  if(hist) raw <-TRUE 
  if(raw) split <- TRUE
  NFACTORS <- nfactors
 
  #check to see if the first parameter is a list of keys, if not, then we want to do reliability on just one scale
  
  if(!is.null(keys)){ if(NCOL(keys)==1) { n.scales <- length(keys) } else {
         items <- keys  #this is the case where the user did not specify keys but just the data
         message("keys not specified, all items will be scored")
         n.scales  <- 1 
         keys <- list()
          keys[["All_items"]] <- colnames(items)
         }} else {n.scales <- 1
   keys[["All_items"]] <- colnames(items)
     }
      if(isCorrelation(items)) {cors<- TRUE
           n.obs=n.obs} else {cors<- FALSE
           n.obs<- NROW(items)}
  #n.obs are needed for CFI in unidim
  for (scales in 1:n.scales) {
  scale.key <- keys[[scales]]
 select <- selectFromKeys(scale.key)
 if(length(select)>1 )  {
  if (length(select)< 4) {nfactors <- 1} else {nfactors <- NFACTORS}   #omegah with 2 factors is not defined for nvar<=5
  if(cors) {om <- omegah(items[select,select], nfactors=nfactors,plot=plot,two.ok=TRUE,fm=fm) } else {
   om <- omegah(items[,select], nfactors=nfactors,n.obs=n.obs,plot=plot,two.ok=TRUE,covar=covar,fm=fm)}
 	uni <- unidim(om$R, n.obs=n.obs)  #use the R object rather than redoing the factoring

 
  if(split){temp.keys <- colnames(om$R) <-  gsub("-","",colnames(om$R))
            sign.key <- rep("",length(select))
            sign.key[which(colnames(om$R)!= rownames(om$R))] <- "-"
             temp.keys <- paste0(sign.key,temp.keys)
                       
    # split.half <- suppressWarnings(splitHalf(om$R,raw=raw,brute=FALSE,n.sample=n.sample, key=temp.keys))
    split.half <- suppressWarnings(splitHalf(om$R,raw=raw,brute=brute,n.sample=n.sample,check.keys=check.keys,covar=covar)) #don't use temp.keys
          best[[scales]] <- list(max=split.half$maxAB)
          worst[[scales]] <- list(min=split.half$minAB)
         
          result[[scales]] <- list(omega_h = om$omega_h, alpha = split.half$alpha, omega.tot = om$omega.tot,
                    u=uni$u[1],tau=uni$u[2],cong=uni$u[3],
                maxrb=split.half$maxrb,minrb=split.half$minrb,
                mean.r=split.half$av.r, med.r <- split.half$med.r, n.items=length(select),
                 CFI <- uni$uni[7], ECV = uni$uni[8])  
          if(raw) {splits[[scales]] <- split.half$raw}} else {      
   result[[scales]] <- list(omega_h = om$omega_h, alpha = om$alpha, omega.tot = om$omega.tot,u=uni$u[1],tau=uni$u[2],
                     cong=uni$u[3],n.items=length(select),
                     CFI <- uni$uni[7], ECV = uni$uni[8])}
  res.name[scales] <- names(keys)[scales]
   }
  }


 best <- unlist(unlist(best,recursive=FALSE),recursive=FALSE) #Creates a list that keeps the names
 worst <- unlist(unlist(worst,recursive=FALSE),recursive=FALSE)
  # names(result) <- res.name
  if(split) {ncol <- 13} else {ncol <- 11}
   result.df <- matrix(unlist(result[!is.null(result)]), ncol=ncol,byrow=TRUE)
  if(split) {   colnames(result.df) <- c("omega_h", "alpha", "omega.tot", "Uni","tau","cong","max.split","min.split","mean.r", "med.r", "n.items","CFI","ECV")} else {
  colnames(result.df) <- c("omega.h", "alpha", "omega.tot","Uni","tau","cong","n.items","CFI","ECV")}
  rownames(result.df) <- unlist(res.name)
  if(raw) {
  lx <- unlist(lapply(splits,length))
  splits <- splits[lx>0]
  names(splits) <- rownames(result.df)
  
 # splits.mat <- matrix(unlist(splits),ncol=length(keys))
 # colnames(splits.mat) <- names(keys)
 
 
   class(result.df) <- c("psych","reliability", "matrix")
   names(best)<- paste(rep(rownames(result.df),each=2),names(best))
   names(worst) <- paste(rep(rownames(result.df),each=2),names(worst))
  result <- list(result.df = result.df,splits= splits,max=best,min=worst, Call = cl)
  if(hist) {multi.hist(splits)}
  class(result) <-  c("psych","reliability")
  return(result) } else {

   class(result.df) <- c("psych","reliability", "matrix")
  return(result.df)
  }
  }
  
  #Created June 11-17, 2021
  #fixed 6/20/21 to avoid the problem of null cases
  #fixed 1/3/23 to avoid problem with single items
  
  
                                                                