

makerepeated <- function(f=NULL, key="STATE",filter=NULL) {
#go get the files we want
if(is.null(f)) {
 fn <- file.choose()  } else { fn <- f}  #find a file in the directory you want
  dir <- dirname(fn)  #the directory where the file was found
  files.list <- list.files(dir)
 
select <- grep(key,files.list,ignore.case=TRUE) #these are the ones that match key
selected <- files.list[select]
if(!is.null(filter)) {selected <- files.list[select]
select <- grep(filter,selected,ignore.case=TRUE,invert=TRUE)
selected <- selected[select]}
#first figure out the n.var for this set of files
for (i in 1:length(selected)) {
file <- selected[i]
temp <- strsplit(file,'[- ]')[[1]]  #strip out the experiment names from the other text
study <- temp[1] #the study name
time <- grep('[1234567890]',temp)
#now, it is possible that the key name is embedded in the time
#time <- gsub(key,"",temp[time],ignore.case=TRUE)
#time <- gsub('[.]',"",temp[time],ignore.case=TRUE)
if(length(time)>0 ){ time <- temp[time] } else {time <- 1}

#we use scan rather than read.fwf because we have different widths for each file and we don't the format
path <- file.path(dir,file)
#this functionally reads fwf and converts to numeric
cat("\nFile read=",file)
temp <- scan(path,what="raw")
periods <- grep("[.]",temp) 
if(length(periods) > 0) {temp <- strsplit(temp,"[.]") #we need to strip out periods
temp <- unlist(temp)} 
#how many records per subject (some files have id field  followed by a code field, followed by the data
if(nchar(temp[2]) < nchar(temp[3])) {nlines  <- 3} else {nlines <-2}  
n.obs <- length(temp)/nlines
n.var <- nchar(temp[nlines])
new.data <- matrix(NA,nrow=n.obs,ncol=n.var +1)

for( j in 1:n.obs) {  #now, form a matrix of characters with the data
new.data[j,1] <- temp[1+(j-1)*nlines]
temp.split <- strsplit(temp[j*nlines],"") 
new.data[j,2:(n.var+1)] <- temp.split[[1]]   #this gets the id field 
}
new.data <- matrix(as.numeric(new.data),ncol=n.var + 1)
colnames(new.data) <- c("id",paste0("V",1:n.var) ) 
new.data.df <- data.frame(study=study,time=time,new.data)  
if(i ==1) {big.data <- new.data.df} else {
if(NCOL(new.data.df) != NCOL(big.data)) {cat("\nOoops  file sizes don't match", file, "\n")} else {

 big.data <- rbind(big.data,new.data.df)}

}
}   #end of for i loop
return(big.data) 
}
