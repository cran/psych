#added  just correlate with criteria to speed it up (June 23, 2017)
"bestItems" <- 
function(x,criteria=1,cut=.1, n.item=10,raw=TRUE, abs=TRUE, dictionary=NULL,check=FALSE,digits=2) {

if(check) {item.var <- apply(x,2,sd,na.rm=TRUE)  #added check 10/14/17
       bad <- which((item.var <= 0)|is.na(item.var))
       if((length(bad) > 0) ) {
            for (baddy in 1:length(bad)) {message( "Item = ",colnames(x)[bad][baddy], " had no variance and was deleted")}
            x <- x[,-bad] 
             }
             }
 result <- list()
for(i in 1:length(criteria)) {criterion <- criteria[i]
if(raw)  { x <- cor(x,x[,criterion],use="pairwise")
    if(NROW(criterion)> 1)  {x <- cbind(x,criterion)   #to allow for a separate object
    criteron <- "criteria" }
    } #the normal case --convert to correlation if necessary
  
if(abs) {ord <- order(abs(x[,criterion]),decreasing=TRUE)
  value <- x[ord,criterion,drop=FALSE]
  count <- sum(abs(value) > cut,na.rm=TRUE)
  if(!is.null(n.item)) count <- max(count,n.item)
  value <- value[1:count,,drop=FALSE]
  } else {ord <- order(x[,criterion],decreasing=TRUE)
  value <- x[ord,criterion]
  value <- value[value,criterion > cut] }
value <- round(data.frame(value),digits)
if((!is.null(dictionary)) && !is.factor(dictionary)) {temp <- lookup(rownames(value),dictionary)
   value <- merge(value,temp,by="row.names",all.x=TRUE,sort=FALSE)
   rownames(value) <- value[,"Row.names"]
   value <- value[,-1]
  if(abs) {ord <- order(abs(value[,criterion]),decreasing=TRUE) } else {ord <- order(value[,criterion],decreasing=TRUE)}
   value <- value[ord,] 
   
   }
result[[criterion]] <- value }
return(result)
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
function(x,y,criteria=NULL) {
if (is.null(criteria)) {temp <- match(x,rownames(y))} else {
     temp <- match(x,y[,criteria])}
     if(any(!is.na(temp))) {
 y <- (y[temp[!is.na(temp)],,drop=FALSE]) } else {y <- NA}
return(y)}
