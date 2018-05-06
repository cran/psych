#Various useful utility functions

 # list the files in a directory holding a particular file, or a particular directory
"filesList" <- function(f=NULL) {
 if(is.null(f)) { f <- file.choose()} else {
   if(dir.exists(f)) {dir <- f } else {dir <- dirname(f)}}  #find a file in the directory you want
  files.list <- list.files(dir)
  cat("\nFiles in the directory", dir, "\n")
  return(files.list)
  }
  
  "filesInfo" <- function(f=NULL,max=NULL) {
   if(is.null(f)) { f <- file.choose()} else {
   if(dir.exists(f)) {dir <- f } else {dir <- dirname(f)}}
    files.list <- list.files(dir)
   if(is.null(max)) max <- length(files.list)
   info <- list(max)
   for(i in 1:max) {
   info[[i]] <- file.info(file.path(dir,files.list[i]))}
   info.df <- info[[1]]
   for (i in 2:max) {
   info.df <- rbind(info.df,info[[i]])}
   info.df <-cbind(file=1:max,info.df)
   return(info.df)
  }
  
  
  "fileScan" <- function(f=NULL,nlines=3,max=NULL,from=1,filter=NULL) {
  if(is.null(f)) {f <- file.choose()}  #find a file in the directory you want
   dir <- dirname(f)  #the directory where the file was found
   files.list <- list.files(dir)
   if(!is.null(filter)) {select <- grep(filter,files.list,ignore.case=TRUE) #these are the ones that match filter            
               files.list <- files.list[select]}
   n.files <- length(files.list)
   if(!is.null(max)) n.files <- max + from
   for (i in from:n.files) {
   file <- files.list[i]
     path <- file.path(dir,file)
     suffix <- file_ext(file)
      if(suffix %in% c("xls","xlsx","doc","sav","data","dat","rds","R","r","RDS", "XPT","xpt","Rda","rda","Rdata","RData","rdata","SYD","syd","sys","jmp","sas7bdat")) {
      cat("\nFile = ",i, "Name = ", file, "Was skipped") } else {
    # temp <- scan(path,what="raw",nlines=nlines)
    temp <- readLines(path,n=nlines)
    cat("\nFile = ",i, "Name = ", file, "\n",temp,"\n")}
    }
    return(dir)
    }
   
  
  #a work around the failure of file.choose(new=TRUE) to work in Rstudio
  
"fileCreate" <- function(newName="new.file") {
  fn <- file.choose()
  dir <- dirname(fn)
  new.path <- file.path(dir,newName) 
  if(!file.exists(new.path)) {
  file.create(new.path)
  return(new.path) } else {cat('\nFile already exists, try a different name')}
  }
  
 #Completely rewritten 1/20/18 to follow  the help pages for order more closely
#sort a data frame according to one or multiple columns 
#will only work for data.frames (not matrices)
dfOrder <- function(object,columns=NULL,absolute=FALSE) {
   if(is.null(columns)) columns <- 1:ncol(object)
	 nc <- length(columns)
	 cn <- colnames(object)
 	 temp <- rep(1,nc)    
 	 if(is.character(columns)) {  #treat character strings 
   		 temp [strtrim(columns,1)=="-"] <- -1
    	 if(any(temp < 0  ) )  {columns <- sub("-","",columns) }
    	 
    	} else {temp[columns < 0] <- -1
    	        columns <- abs(columns) }
    	
   if(is.character(columns) ) {	 for (i in 1:length(columns)) {columns[i] <- (which(colnames(object) == columns[i]))
       }
       }
      
   	  columns <- colnames(object)[as.numeric(columns)]
    if(absolute) {  temp.object<- t(t(abs(char2numeric(object[columns]))) * temp)  } else {
   	  temp.object<- t(t(char2numeric(object[columns])) * temp)}
   	  temp.object <- data.frame(temp.object)
 
   	 ord <- do.call(order,temp.object)
     if(length(ord) > 1) {
   	   return(object[ord,]) }else {return(object)} #added length test 4/26/18
       }

