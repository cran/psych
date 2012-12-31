"df2latex" <- 
function(x,digits=2,rowlabels=TRUE,apa=TRUE,short.names=TRUE,font.size ="scriptsize", heading="A table from R",caption="df2latex") {
#first set up the table
 nvar <- dim(x)[2]
comment <- paste("%", match.call())
header <- paste("\\begin{",font.size,"} \\begin{table}[htdp]",
"\\caption{",caption,"}
\\begin{center}
\\begin{tabular}",sep="")
if(rowlabels) {header <- c(header,"{l",rep("r",(nvar)),"}\n")} else {header <- c(header,"{",rep("r",(nvar)),"}\n")}
if(apa) {header <- c(header,
"\\multicolumn{",nvar,"}{l}{",heading,"}",
'\\cr \n \\hline ')
footer <- paste(" \\hline ")}  else {footer <- NULL}
footer <- paste(footer,"
\\end{tabular}
\\end{center}
\\label{default}
\\end{table} 
\\end{",font.size,"}
",sep=""
)

#now put the data into it
if(!is.null(digits)) {if(is.numeric(x) ) {x <- round(x,digits=digits)} else {x <- try(round(x,digits=digits)) }}
 
 cname <- colnames(x)
 if (short.names) cname <- 1:nvar
 names1 <- paste(cname[1:(nvar-1)], " & ")
 lastname <- paste(cname[nvar],"\\cr \n")
 
 if(apa)  {allnames <- c("Variable  &  ",names1,lastname," \\hline \n")} else {if(rowlabels) {allnames <- c("  &  ",names1,lastname,"\\cr \n")} else {
             allnames <- c(names1,lastname,"\\cr \n")}}
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
 
 
"cor2latex" <- function(x,digits=2,rowlabels=TRUE,lower=TRUE,apa=TRUE,short.names=TRUE,font.size ="scriptsize",heading="A correlation table from R",caption="cor2latex") {
if(nrow(x) > ncol(x) ) {x <- cor(x,use="pairwise")}
r <- round(x,digits)
r <- format(r,nsmall=digits)  #this converts to character but keeps the right number of digits
if(lower) {r[upper.tri(r)] <- "~"} else {r[lower.tri(r)] <- "~"} 
return(df2latex(r,digits=NULL,rowlabels=rowlabels,apa=apa,short.names=short.names,font.size=font.size,heading=heading,caption=caption))
}

"fa2latex" <- 
function(f,digits=2,rowlabels=TRUE,apa=TRUE,short.names=FALSE,font.size ="tiny", heading="A factor analysis table from R",caption="fa2latex") {
x <- unclass(f$loadings)
if(!is.null(f$Phi)) {Phi <- f$Phi} else {Phi <- NULL}
nfactors <-ncol(x)

if(nfactors > 1) {if(is.null(Phi)) {h2 <- rowSums(x^2)} else {h2 <- diag(x %*% Phi %*% t(x)) }} else {h2 <-x^2}
u2 <- 1- h2
vtotal <- sum(h2 + u2)
x <- data.frame(x,h2=h2,u2=u2)
#first set up the table
 nvar <- dim(x)[2]
comment <- paste("%", match.call())
header <- paste("\\begin{",font.size,"} \\begin{table}[htdp]",
"\\caption{",caption,"}
\\begin{center}
\\begin{tabular}",sep="")
header <- c(header,"{l",rep("r",nvar),"}\n")
if(apa) header <- c(header,
"\\multicolumn{",nvar,"}{l}{",heading,"}",
'\\cr \n \\hline ')
if(apa) {footer <- paste(" \\hline ")} 
footer <- paste(footer,"
\\end{tabular}
\\end{center}
\\label{default}
\\end{table} 
\\end{",font.size,"}
",sep=""
)


#now put the data into it
 x <- round(x,digits=digits)
 cname <- colnames(x)
 if (short.names) cname <- 1:nvar
 names1 <- paste(cname[1:(nvar-1)], " & ")
 lastname <- paste(cname[nvar],"\\cr \n")
 
 if(apa)  {allnames <- c("Variable  &  ",names1,lastname," \\hline \n")} else {allnames <- c("  &  ",names1,lastname,"\\cr \n")}
 x <- format(x)  #to keep the digits the same
 value <- apply(x,1,paste,collapse="  &  ") #insert & between columns
 if(rowlabels) value <- paste(names(value),"  & ",value)
 values <- paste(value, "\\cr", "\n")  #add \\cr at the end of each row

 #now put it all together
 cat(comment,"\n")  #a comment field saying where the data came from
 cat(header)   #the header information
 cat(allnames) #the variable names
 cat(values)  #the factor loadings
 
 
 #now find and show the variance accounted for
 x <- f$loadings     #use the original values
 nvar <- nrow(x)
  if(is.null(Phi)) {if(nfactors > 1)  {vx <- colSums(x^2) } else {
                                      vx <- diag(t(x) %*% x)
                                      vx <- vx*nvar/vtotal 
      	                             }} else {vx <- diag(Phi %*% t(x) %*% x)
      	                                      vx <- vx*nvar/vtotal }
      	  #names(vx) <- colnames(x)[1:nvar]
      	  vx <- round(vx,digits) 
          loads <- c("\\hline \\cr SS loadings &",paste(vx," & ",sep=""))
          cat(loads)
           
          #varex <- rbind("SS loadings " =   vx)
         # varex <- rbind(varex, "Proportion Var" =  vx/nvar)
          # if (nfactors > 1) {varex <- rbind(varex, "Cumulative Var"=  cumsum(vx/nvar))
          #varex <- rbind(varex, "Cum. factor Var"=     cumsum(vx/sum(vx)))}
    
 if(!is.null(Phi)) {
        cat("\\cr 
            \\hline \\cr \n")
        Phi <- round(Phi,digits)
        phi <- format(Phi,nsmall=digits)
       phi <-apply(phi,1,paste,collapse=" & ")
       phi <-paste(colnames(x),"  &",phi)
       phi <- paste(phi, "\\cr", "\n")
       cat(phi)}
 cat(footer)   #close it up with a footer

 }
 
 
 
 "irt2latex" <- 
function(f,digits=2,rowlabels=TRUE,apa=TRUE,short.names=FALSE,font.size ="tiny", heading="An IRT factor analysis table from R",caption="fa2latex") {
i <- 1  
x <- f$plot$sumInfo[[i]]

#first set up the table
 nvar <- ncol(x)
comment <- paste("%", match.call())
header <- paste("\\begin{",font.size,"} \\begin{table}[htdp]",
"\\caption{",caption,"}
\\begin{center}
\\begin{tabular}",sep="")
header <- c(header,"{l",rep("r",nvar),"}\n")
if(apa) header <- c(header,
"\\multicolumn{",nvar,"}{l}{",heading,"}",
"\\cr  \\hline \\cr",
"\n & \\multicolumn{7}{c}{Item information at $\\theta$}  \\cr \\cline{2-8}  ")
if(apa) {footer <- paste(" \\hline ")} 
footer <- paste(footer,"
\\end{tabular}
\\end{center}
\\label{default}
\\end{table} 
\\end{",font.size,"}
",sep=""
)


#now put the data into it
 x <- round(x,digits=digits)
 cname <- colnames(x)
 if (short.names) cname <- 1:nvar
 names1 <- paste(cname[1:(nvar-1)], " & ")
 lastname <- paste(cname[nvar],"\\cr \n")
 
 if(apa)  {allnames <- c("Item  &  ",names1,lastname," \\hline \n")} else {allnames <- c("  &  ",names1,lastname,"\\cr \n")}
 x <- format(x)  #to keep the digits the same
 value <- apply(x,1,paste,collapse="  &  ") #insert & between columns
 if(rowlabels) value <- paste(names(value),"  & ",value)
 values <- paste(value, "\\cr", "\n")  #add \\cr at the end of each row



 #now put it all together
 cat(comment,"\n")  #a comment field saying where the data came from
 cat(header)   #the header information
 cat(allnames) #the variable names
 cat(values)  #the item information 
 cat("\\hline \n & \\multicolumn{7}{c}{Summary statistics at $\\theta$} \\cr \\cline{2-8}")
 test.info <- colSums(f$plot$sumInfo[[i]])
 sem <- sqrt(1/test.info)
 reliab <- 1 - 1/test.info
 summary <- rbind(test.info,sem,reliab)
 summary <- round(summary,digits)
 summary <- format(summary,nsmall=digits)
 summary <- cbind(c("Test.info","SEM","Reliability"),summary)
 summary <- apply(summary,1,paste,collapse="  & ")
 summary <- paste(summary,"\\cr \n") 
 cat(summary)
 
 cat(footer)   #close it up with a footer

 }
 
 

