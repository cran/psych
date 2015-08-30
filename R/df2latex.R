#modified April 6, 2015 to return the table invisibly as well so it can be embedded in a Sweave document
#November 22, 2013  Modified with help from Davide Morselli to allow for "stars"
#also allows for printing straight text (char=TRUE)
#cor2latex was modified following Davide Morselli's suggestion to allow direct calculation of the correlations
#
"df2latex" <- 
function(x,digits=2,rowlabels=TRUE,apa=TRUE,short.names=TRUE,
font.size ="scriptsize",big.mark=NULL, drop.na=TRUE, heading="A table from the psych package in R",
caption="df2latex",label="default",char=FALSE,stars=FALSE,silent=FALSE,file=NULL,append=FALSE) {
#first set up the table
 nvar <- dim(x)[2]
 rname<- rownames(x)
comment <- paste("%", match.call())
header <- paste("\\begin{table}[htpb]",
"\\caption{",caption,"}
\\begin{center}
\\begin{",font.size,"} 
\\begin{tabular}",sep="")
if(stars) {if(rowlabels) {
               header <- c(header,"{l",rep("S",(nvar)),"}\n")} else {header <- c(header,"{",rep("S",(nvar+1)),"}\n")}  } else {
              if(rowlabels) { header <- c(header,"{l",rep("r",(nvar)),"}\n")} else {header <- c(header,"{",rep("r",(nvar+1)),"}\n")}
               }
if(apa) {header <- c(header,
"\\multicolumn{",nvar,"}{l}{",heading,"}",
'\\cr \n \\hline ')
footer <- paste(" \\hline ")}  else {footer <- NULL}
if (stars){
      footer <- paste(" \\hline 
                      \n \\multicolumn{7}{l}{\\scriptsize{\\emph{Note: }\\textsuperscript{***}$p<.001$; 
                      \\textsuperscript{**}$p<.01$; 
                      \\textsuperscript{*}$p<.05$",".}}" ,sep = "")
    }else{
      footer <- paste(" \\hline ")}
footer <- paste(footer,"
\\end{tabular}
\\end{",font.size,"}
\\end{center}
\\label{",label,"}
\\end{table} 

",sep=""
)

#now put the data into it
if(!char) {if(!is.null(digits)) {if(is.numeric(x) ) {x <- round(x,digits=digits)} else {x <- try(round(x,digits=digits)) }}
      }
 
 cname <- colnames(x)
 if (short.names) cname <- abbreviate(cname,minlength=digits+3)  #cname <- 1:nvar

 names1 <- paste(cname[1:(nvar-1)], " & ")
 lastname <- paste(cname[nvar],"\\cr \n")
 
if(apa)  {allnames <- c("Variable  &  ",names1,lastname," \\hline \n")} else {if(rowlabels) {allnames <- c("  &  ",names1,lastname,"\\cr \n")} else {
             allnames <- c(names1,lastname,"\\cr \n")}}
if(!char) {if(is.null(big.mark)) { x <- format(x,drop0trailing=FALSE)} else   #to keep the digits the same
                      {x <- prettyNum(x,big.mark=",",drop0trailing=FALSE)} 
   }   
 value <- apply(x,1,paste,collapse="  &  ") #insert & between columns

 if(rowlabels) {value <- paste(sanitize.latex(rname),"  & ",value)} else {value <- paste("  & ",value)}
 values <- paste(value, "\\cr", "\n")  #add \\cr at the end of each row
 if(drop.na) values <- gsub("NA","  ",values,fixed=TRUE)

 #now put it all together
 if(!silent) {cat(comment,"\n")  #a comment field saying where the data came from
 cat(header)   #the header information
 cat(allnames) #the variable names
 cat(values)  #the data
 cat(footer)   #close it up with a footer
 } 
result <- c(header,allnames,values,footer)
if(!is.null(file)) write.table(result,file=file,row.names=FALSE,col.names=FALSE,quote=FALSE,append=append)

invisible(result)
 }
 
 
 
 cor2latex <- function (x, use = "pairwise", method="pearson", adjust="holm", stars = FALSE, digits=2, rowlabels = TRUE, lower = TRUE, apa = TRUE, 
                       short.names = TRUE, font.size = "scriptsize", heading = "A correlation table from the psych package in R.", 
                       caption = "cor2latex", label = "default",silent=FALSE,file=NULL,append=FALSE)
{
if(stars) heading  <- paste(heading, "Adjust for multiple tests = ",adjust )
if (!is.na(class(x)[2]) & class(x)[2]=="corr.test") {  #we already did the analysis, just report it
      r <- x$r
      p <- x$p} else {
 
      if (nrow(x) > ncol(x)) {   #find the correlations 
        x <- corr.test(x, use=use,method=method,adjust=adjust)
        r <- x$r
        p <- x$p   } else {   #take the correlations as given
        r <- x
        p <- NULL
      }
    }
    r <- round(r, digits)
    r <- format(r, nsmall = digits,drop0trailing=FALSE)  #this converts to character but keeps the right number of digits)
    if (lower) {
      r[upper.tri(r)] <- "~"
    } else {
      r[lower.tri(r)] <- "~"
    }
  
    if(stars && is.null(p))  stop("To print significance levels, x must be be either a data frame of observations or a correlation matrix created with the corr.test function of the package psych. If you are not interested in displaying signicance level set stars = FALSE")
     
          #p[upper.tri(p,diag=FALSE)]  #the adjusted probability values
       
 
    mystars <- ifelse(p < .001, "{***}", ifelse(p < .01, "{**}", ifelse(p < .05, "{*}", "")))
    mystars <- t(mystars)
   if(stars) { R <- matrix(paste(r,mystars,sep=""),ncol=ncol(r))} else {R <- r}
    diag(R) <- paste(diag(r), " ", sep="")
    rownames(R) <- colnames(r)
    colnames(R) <- colnames(r)
    if (lower) {
      R[upper.tri(R, diag = FALSE)] <- ""
     
    } else {
      R[lower.tri(R, diag = FALSE)] <- ""
      
    }
     
    if(stars) {char<- TRUE} else {char <- FALSE}
  return(df2latex(R, digits = digits, rowlabels = rowlabels, 
                  apa = apa, short.names = short.names, font.size = font.size, 
                  heading = heading, caption = caption, label = label, char=TRUE,stars = stars,silent=silent,file=file,append=append))
}

"fa2latex" <- 
function(f,digits=2,rowlabels=TRUE,apa=TRUE,short.names=FALSE,cumvar=FALSE,cut=0,big=.3,alpha=.05,font.size ="scriptsize", heading="A factor analysis table from the psych package in R",caption="fa2latex",label="default",silent=FALSE,file=NULL,append=FALSE) {
if(class(f)[2] == "fa.ci") {
if(is.null(f$cip)) {px <- f$cis$p} else {px <- f$cip}} else {px <- NULL}  #get the probabilities if we did fa.ci
#if(class(f)[2] !="fa") f <- f$fa
x <- unclass(f$loadings)
if(!is.null(f$Phi)) {Phi <- f$Phi} else {Phi <- NULL}
nfactors <- ncol(x)

if(nfactors > 1) {if(is.null(Phi)) {h2 <- rowSums(x^2)} else {h2 <- diag(x %*% Phi %*% t(x)) }} else {h2 <-x^2}
u2 <- 1- h2
vtotal <- sum(h2 + u2)
if(cut > 0) x[abs(x) < cut] <- NA    #modified May 13 following a suggestion from Daniel Zingaro
if(!is.null(f$complexity)) {x <- data.frame(x,h2=h2,u2=u2,com=f$complexity) } else {x <- data.frame(x,h2=h2,u2=u2)}
#first set up the table
 nvar <- dim(x)[2]
comment <- paste("% Called in the psych package ", match.call())
header <- paste("\\begin{table}[htpb]",
"\\caption{",caption,"}
\\begin{center}
\\begin{",font.size,"} 
\\begin{tabular}",sep="")
header <- c(header,"{l",rep("r",nvar),"}\n")
if(apa) header <- c(header,
"\\multicolumn{",nvar,"}{l}{",heading,"}",
'\\cr \n \\hline ')
if(apa) {footer <- paste(" \\hline ")} 
footer <- paste(footer,"
\\end{tabular}
\\end{",font.size,"}
\\end{center}
\\label{",label,"}
\\end{table} 

",sep=""
)


#now put the data into it
 
 x <- round(x,digits=digits)   
 
 
 cname <- colnames(x)
 if (short.names) cname <- 1:nvar
 names1 <- paste(cname[1:(nvar-1)], " & ")
 lastname <- paste(cname[nvar],"\\cr \n")
 
 if(apa)  {allnames <- c("Variable  &  ",names1,lastname," \\hline \n")} else {allnames <- c("  &  ",names1,lastname,"\\cr \n")}
 fx <- format(x,drop0trailing=FALSE)  #to keep the digits the same
 {if(!is.null(px) && (cut == 0)) { temp <- fx[1:nfactors]
 temp[px < alpha] <- paste("\\bf{",temp[px < alpha],"}",sep="")
 fx[1:nfactors] <- temp
 }
 if(big > 0) {temp <- fx[1:nfactors]  
   x <- x[1:nfactors]
  temp[!is.na(x) & (abs(x) > big)] <- paste("\\bf{",temp[!is.na(x) & (abs(x) > big)],"}",sep="")
   fx[1:nfactors] <- temp
   }
 value <- apply(fx,1,paste,collapse="  &  ") #insert & between columns
 value <- gsub("NA", "  ", value, fixed = TRUE)
 if(rowlabels) value <- {paste(sanitize.latex(names(value)),"  & ",value)} else {paste("  &  ",value)}
 values <- paste(value, "\\cr", "\n")  #add \\cr at the end of each row

 #now put it all together
 if(!silent) {
 cat(comment,"\n")  #a comment field saying where the data came from
 cat(header)   #the header information
 cat(allnames) #the variable names
 cat(values)  #the factor loadings
 }
 
 #now find and show the variance accounted for
 x <- f$loadings     #use the original values not the rounded ones
 nvar <- nrow(x)
  if(is.null(Phi)) {if(nfactors > 1)  {vx <- colSums(x^2) } else {
                                      vx <- diag(t(x) %*% x)
                                      vx <- vx*nvar/vtotal 
      	                             }} else {vx <- diag(Phi %*% t(x) %*% x)
      	                                      vx <- vx*nvar/vtotal }
      	  #names(vx) <- colnames(x)[1:nvar]
      	  vx <- round(vx,digits) 
      	  
          loads <- c("\\hline \\cr SS loadings &",paste(vx," & ",sep=""),"\\cr  \n")
           
 if(!silent) { cat(loads)}
       summ <- NULL
            
          #varex <- rbind("SS loadings " =   vx)
          if(cumvar) {
          provar <- round(vx/nvar,digits)        
         summ <- c("Proportion Var &" ,paste(  provar, "  & ",sep=""),"\\cr \n")
         
       #  cat("Proportion Var &" ,paste(  provar, "  & ",sep=""),"\\cr \n")
           if (nfactors > 1) {cumvar <- round(cumsum(vx/nvar),digits)
             cumfavar <- round(cumsum(vx/sum(vx)),digits=digits)
        summ <- c(summ,  "Cumulative Var & ",paste( cumvar," & ", sep=""),"\\cr \n",
         "Cum. factor Var & ",paste(round(cumsum(vx/sum(vx)),digits=digits),"  & ",sep=""),"\\cr \n")
          } 
          if(!silent) {cat(summ)  }
          }
   loads <- c(loads,summ)      
 if(!is.null(Phi)) {
      summ <-   c("\\cr 
            \\hline \\cr \n") 
            if(!silent) {cat(summ)  }
        Phi <- round(Phi,digits)
        phi <- format(Phi,nsmall=digits)
       phi <-apply(phi,1,paste,collapse=" & ")
       phi <-paste(colnames(x),"  &",phi)
       phi <- paste(phi, "\\cr", "\n")
       loads <- c(loads,summ,phi)
     if(!silent) {  cat(phi)}
     }
if(!silent) { cat(footer)}   #close it up with a footer
 }
 values <- c(values,loads)
 result <- c(header,allnames,values,footer)
 if(!is.null(file)) write.table(result,file=file,row.names=FALSE,col.names=FALSE,quote=FALSE,append=append)

invisible(result)
 }
 
 
 "irt2latex" <- 
function(f,digits=2,rowlabels=TRUE,apa=TRUE,short.names=FALSE,font.size ="scriptsize", heading="An IRT factor analysis table from R",caption="fa2latex" ,label="default",silent=FALSE,file=NULL,append=FALSE) {
nf <- length(f$plot$sumInfo)   #create nf tables 
for(i in (1:nf)) { 
x <- f$plot$sumInfo[[i]]
if(nf>1) {
rowmax <- apply(x,1,max, na.rm=TRUE)
 rowmax <- which(rowmax <.001,arr.ind=TRUE) 
 if(!is.null(rowmax)) x <- x[-rowmax,]}
#first set up the table
 nvar <- ncol(x)
comment <- paste("%", match.call())
header <- paste("\\begin{",font.size,"} \\begin{table}[htpb]",
"\\caption{",caption,"}
\\begin{center}
\\begin{tabular}",sep="")
header <- c(header,"{l",rep("r",nvar),"}\n")
if(apa) header <- c(header,
"\\multicolumn{",nvar,"}{l}{",heading," for factor " , i, " }",
"\\cr  \\hline \\cr",
"\n & \\multicolumn{7}{c}{Item information at $\\theta$}  \\cr \\cline{2-8}  ")
if(apa) {footer <- paste(" \\hline ")} 
footer <- paste(footer,"
\\end{tabular}
\\end{center}
\\label{",label,"}
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
 x <- format(x,drop0trailing=FALSE)  #to keep the digits the same
 value <- apply(x,1,paste,collapse="  &  ") #insert & between columns
 if(rowlabels) value <- paste(sanitize.latex(names(value)),"  & ",value)
 values <- paste(value, "\\cr", "\n")  #add \\cr at the end of each row

 #now put it all together

 test.info <- colSums(f$plot$sumInfo[[i]])
 sem <- sqrt(1/test.info)
 reliab <- 1 - 1/test.info
 summary <- rbind(test.info,sem,reliab)
 summary <- round(summary,digits)
 summary <- format(summary,nsmall=digits)
 summary <- cbind(c("Test.info","SEM","Reliability"),summary)
 summary <- apply(summary,1,paste,collapse="  & ")
 summary <- paste(summary,"\\cr \n") 
 
if(!silent) {  cat(comment,"\n")  #a comment field saying where the data came from
 cat(header)   #the header information
 
 cat(allnames) #the variable names
 cat(values)  #the item information 
 cat("\\hline \n & \\multicolumn{7}{c}{Summary statistics at $\\theta$} \\cr \\cline{2-8}")
 cat(summary)
 cat(footer)   #close it up with a footer'
 }
 }
 result <- c(header,allnames,values,summary,footer)
 if(!is.null(file)) write.table(result,file=file,row.names=FALSE,col.names=FALSE,quote=FALSE,append=append)

invisible(result)
 }
 
 
 #adapted from various sources, including xtable
"sanitize.latex" <- 
function(astring) {
result <- astring
result <- gsub("&", "\\&", result, fixed = TRUE)
result <- gsub("_", "\\_", result, fixed = TRUE)
result <- gsub("%", "\\%", result, fixed = TRUE)
return(result)
}
 

 
 
 #added December 28, 2013
 "omega2latex" <- 
function(f,digits=2,rowlabels=TRUE,apa=TRUE,short.names=FALSE,cumvar=FALSE,cut=.2,font.size ="scriptsize", heading="An omega analysis table from the psych package in R",caption="omega2latex",label="default",silent=FALSE,file=NULL,append=FALSE) {
if(class(f)[2] == "omega" ) f$loadings <- f$schmid$sl
x <- unclass(f$loadings)

nfactors <- ncol(x)

h2 <- rowSums(x^2)
u2 <- 1- h2
vtotal <- sum(h2 + u2)

#first set up the table
 nvar <- dim(x)[2]
comment <- paste("% Called in the psych package ", match.call())
header <- paste("\\begin{",font.size,"} \\begin{table}[htpb]",
"\\caption{",caption," with cut = ",cut,"\n $\\omega_h  = ",round(f$omega_h,digits), "\\;\\;\\;\\alpha (\\lambda_3) = ",round(f$alpha,digits), "\\;\\;\\;\\lambda_6^* = ",round(f$G6,digits),"\\;\\;\\; \\omega_t = ",round(f$omega.tot,digits),"$ }
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
\\label{",label,"}
\\end{table} 
\\end{",font.size,"}
",sep=""
)

#now put the data into it
 x[abs(x) < cut] <- NA
 x <- round(x,digits=digits)
 cname <- colnames(x)
 if (short.names) cname <- 1:nvar
 names1 <- paste(cname[1:(nvar-1)], " & ")
 lastname <- paste(cname[nvar],"\\cr \n")
 
 if(apa)  {allnames <- c("Variable  &  ",names1,lastname," \\hline \n")} else {allnames <- c("  &  ",names1,lastname,"\\cr \n")}
 x <- format(x,drop0trailing=FALSE)  #to keep the digits the same
 
  
 value <- apply(x,1,paste,collapse="  &  ") #insert & between columns
 value <- gsub("NA", "  ", value, fixed = TRUE)
 if(rowlabels) value <- {paste(sanitize.latex(names(value)),"  & ",value)} else {paste("  &  ",value)}
 values <- paste(value, "\\cr", "\n")  #add \\cr at the end of each row

 #now put it all together

 
 #now find and show the variance accounted for
 x <- f$loadings     #use the original values
 nvar <- nrow(x)
vx <- colSums(x^2)[1:(ncol(x)-3)]
vx <- round(vx,digits) 
loads <- c("\\hline \\cr SS loadings &",paste(vx," & ",sep=""),"\\cr  \n")

if(!silent) { cat(comment,"\n")  #a comment field saying where the data came from
 cat(header)   #the header information
 cat(allnames) #the variable names
 cat(values)  #the factor loadings
cat(loads)
cat(footer)   #close it up with a footer
}
result <- c(header,allnames,values,loads,footer)
if(!is.null(file)) write.table(result,file=file,row.names=FALSE,col.names=FALSE,quote=FALSE,append=append)

invisible(result)
 }



#added 1/6/14
"ICC2latex" <- 
function(icc,digits=2,rowlabels=TRUE,apa=TRUE,ci=TRUE,
font.size ="scriptsize",big.mark=NULL, drop.na=TRUE, heading="A table from the psych package in R",
caption="ICC2latex",label="default",char=FALSE,silent=FALSE,file=NULL,append=FALSE) {
if((length(class(icc)) < 2 ) | (class(icc)[2] !="ICC")) icc <- ICC(icc)  #do the analysis in case we have not done it yet
#first set up the table
x <- icc$results
 nvar <- dim(x)[2]
 rname<- rownames(x)
comment <- paste("%", match.call())
header <- paste("\\begin{",font.size,"} \\begin{table}[[htpb]",
"\\caption{",caption,"}
\\begin{tabular}",sep="")

if(rowlabels) { header <- c(header,"{l",rep("r",(nvar)),"}\n")} else {header <- c(header,"{",rep("r",(nvar+1)),"}\n")
               }
if(apa) {header <- c(header,
"\\multicolumn{",5,"}{l}{",heading,"}", '\\cr \n \\hline ')
footer <- paste(" \\hline \\cr \\multicolumn{ 5 }{c}{   Number of subjects = ", icc$n.obs, "Number of raters = ",icc$n.judge,"}")}  else {footer <- NULL}


footer <- paste(footer,"
\\end{tabular}
\\label{",label,"}
\\end{table} 
\\end{",font.size,"}
",sep=""
)

#now put the data into it

x[2:nvar] <- try(round(x[2:nvar],digits=digits)) 
    
 
 cname <- colnames(x)
 if(!ci) nvar <- nvar-2
 names1 <- paste(cname[1:(nvar-1)], " & ")
 lastname <- paste(cname[nvar],"\\cr \n")
 
if(apa)  {allnames <- c("Variable  &  ",names1,lastname," \\hline \n")} else {if(rowlabels) {allnames <- c("  &  ",names1,lastname,"\\cr \n")} else {
             allnames <- c(names1,lastname,"\\cr \n")}}
if(!char) {if(is.null(big.mark)) { x <- format(x[1:nvar],drop0trailing=FALSE)} else   #to keep the digits the same
                      {x <- prettyNum(x,big.mark=",",drop0trailing=FALSE)} 
   }   
 value <- apply(x,1,paste,collapse="  &  ") #insert & between columns

 if(rowlabels) {value <- paste(sanitize.latex(rname),"  & ",value)} else {value <- paste("  & ",value)}
 values <- paste(value, "\\cr", "\n")  #add \\cr at the end of each row
 if(drop.na) values <- gsub("NA","  ",values,fixed=TRUE)

 #now put it all together
if(!silent) { cat(comment,"\n")  #a comment field saying where the data came from
 cat(header)   #the header information
 cat(allnames) #the variable names
 cat(values)  #the data
 cat(footer)   #close it up with a footer
 }
 result <- c(header,allnames,values,footer)
 if(!is.null(file)) write.table(result,file=file,row.names=FALSE,col.names=FALSE,quote=FALSE,append=append)

invisible(result)

 }
 

