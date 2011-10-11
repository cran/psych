"df2latex" <- 
function(x,digits=2,rowlabels=TRUE,caption="df2latex") {
#first set up the table
 nvar <- dim(x)[2]
comment <- paste("%", match.call())
header <- c( "\\begin{table}[htdp]",
"\\caption{",caption,"}
\\begin{center}
\\begin{tabular}")
header <- c(header,"{",rep("r",(nvar+ 1)),"}\n")
footer <- paste("\\end{tabular}
\\end{center}
\\label{default}
\\end{table}%")

#now put the data into it
 if(is.numeric(digits) ) x <- round(x,digits=digits)
 names1 <- paste(colnames(x)[1:(nvar-1)], " & ")
 lastname <- paste(colnames(x)[nvar],"\\cr \n")
 allnames <- c("  &  ",names1,lastname,"\\cr \n")
 x <- format(x)  #to keep the digits the same
 value <- apply(x,1,paste,collapse="  &  ") #insert & between columns
 if(rowlabels) value <- paste(names(value),"  & ",value)
 values <- paste(value, "\\cr", "\n")  #add \\cr at the end of each row

 #now put it all together
 cat(comment,"\n")  #a comment field saying where the data came from
 cat(header)   #the header information
 cat(allnames) #the variable names
 cat(values)  #the data
 cat(footer)   #close it up with a footer
 }