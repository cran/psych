#added  just correlate with criteria to speed it up (June 23, 2017)
#fixed 7/16/22 to correctly handle the problem of raw data
"bestItems" <- 
function(x,criteria=1,cut=.1, n.item=10, abs=TRUE, dictionary=NULL,check=FALSE,digits=2,use="pairwise",method="pearson") {

if(check) {item.var <- apply(x,2,sd,na.rm=TRUE)  #added check 10/14/17
       bad <- which((item.var <= 0)|is.na(item.var))
       if((length(bad) > 0) ) {
            for (baddy in 1:length(bad)) {message( "Item = ",colnames(x)[bad][baddy], " had no variance and was deleted")}
            x <- x[,-bad] 
             }
             }
 result <- list()
 best <- list()
 key <- list()
 if(isCorrelation(x)) {raw=FALSE
   external<- FALSE
   cn.crit <- criteria} else 
  {data <- x   #find the correlations 
   if(NROW(criteria)==NROW(data)) {cn.crit <- colnames(criteria)
     data <- cbind(data,criteria)
     external <- TRUE
     criteria <- cn.crit} else {cn.crit <- criteria
      external <- FALSE}
    x <- cor(data,data[,criteria,drop=FALSE],use=use,method=method)
    }  
  if(external) {x[criteria,criteria] <- 0 } 
 #next part removed 7/17/22    
for(i in 1:length(criteria)) {criterion <- criteria[i]
# if(raw)  { x <- cor(data,data[,criterion,drop=FALSE],use="pairwise")
#     if(NROW(criterion)> 1)  {x <- cbind(x,criterion)   #to allow for a separate object
#     criteron <- "criteria" }
#     colnames(x)<- criterion
#     } #the normal case --convert to correlation if necessary

x[criterion,criterion] <- 0
if(abs) {ord <- order(abs(x[,criterion,drop=FALSE]),decreasing=TRUE)
value <- x[ord,criterion,drop=FALSE]
  count <- sum(abs(value) > cut,na.rm=TRUE)
  if(!is.null(n.item)) count <- min(count,n.item)
  value <- value[1:count,,drop=FALSE]
  } else {ord <- order(x[,criterion,drop=FALSE],decreasing=TRUE)
  value <- x[ord,criterion,drop=FALSE]
  value <- value[value,criterion > cut,drop=FALSE] }

  
value <- round(data.frame(value),digits)
colnames(value) <- cn.crit[i]    #this is a kludge to get around a problem with data.frame renaming the variable

if((!is.null(dictionary)) && !is.factor(dictionary)) {temp <- lookup(rownames(value),dictionary)
   value <- merge(value,temp,by="row.names",all.x=TRUE,sort=FALSE)
   rownames(value) <- value[,"Row.names"]
   value <- value[,-1]
  
  if(abs) {ord <- order(abs(value[,criterion]),decreasing=TRUE) } else {ord <- order(value[,criterion],decreasing=TRUE)}
   value <- value[ord,] 
   }

best[[criterion]] <- value }
keys <- best2list(best)
result <- list(bestItems=best,bestKeys=keys)
return(result)
}

best2list  <-function(values) {
keys.list <- list()
nkeys <- length(values)
for (i in 1:nkeys) {
 temp <- rownames(values[[i]] )
 rev <- which (values[[i]][1] < 0)
   if(length(rev)  > 0 ) { temp[rev] <- paste0("-",temp[rev]) }
keys.list[[i]] <- temp
}
names(keys.list) <- names(values)
return(keys.list)
}


#adjusted 11/15/20 to add correlations if provided
 "lookupFromKeys" <- 
function(keys.list,dictionary,n=20,cors = NULL,sort=TRUE,suppress.names=FALSE,digits=2){
n.scales <- length(keys.list)
results <- item.cors <- result.df <- list()
for(i in 1:n.scales) {
  list.name <- names(keys.list[i])
  list.i <- keys.list[[i]]
   keys <- rep(1,length(list.i))[1:(min(n,length(list.i)))]
    neg <- grep("-", list.i[1:(min(n,length(list.i)))])
    keys[neg] <- -1
  select <- sub("-", "", list.i)
  results[[i]] <- lookup(select[1:(min(n,length(list.i)))],dictionary)
  
 if(!is.null(rownames(results[[i]])[keys < 0]))  rownames(results[[i]])[keys < 0] <- paste0(rownames(results[[i]])[keys<0],"-")

 if(!is.null(cors)) { item.cors[[i]] <- round(cors[select[1:(min(n,length(select)))],i],digits=digits)
 result.df[[i]] <- data.frame(results[[i]],cors=item.cors[[i]])
 if(sort) {
      ord <- order(abs(item.cors[[i]]),decreasing=TRUE)
      result.df[[i]] <- result.df[[i]][ord,]
 }
 } else {result.df[[i]] <- data.frame(results[[i]])} #results[[i]] <- c(results[[i]],cors= round(cors[select[1:n],i],digits=digits))
  if(suppress.names) names(results[[i]]) <- ""
 
 # names(results[i]) <- list.name
}
names(result.df) <- names(keys.list)
return(result.df)}  
  
  #lookup which x's are found in y[c1],return matches for y[]
 
 "lookup" <- 
function(x,y,criteria=NULL,keep.na=FALSE) {
if (is.null(criteria)) {temp <- match(x,rownames(y))} else {
     temp <- match(x,y[,criteria])}
     if(any(!is.na(temp))) {	
 y <- (y[temp[!is.na(temp)],,drop=FALSE]) } else {y <- NA}
return(y)}
