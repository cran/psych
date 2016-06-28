#A function to report the difference between two factor models
#adapted from John Fox's sem anova 

anova.psych <- function(object,object2,...) {

name.1 <- deparse(substitute(object))
name.2 <- deparse(substitute(object2))
object1 <- object  #rather cluncky, but gets around the problem of generic anova
#check if we have omega input
if(class(object1)[2] =="omega") object1 <- object1$schmid
if(class(object2)[2] =="omega") object2 <- object2$schmid

chi.1 <- object1$STATISTIC
dof.1  <- object1$dof
BIC.1  <- object1$BIC
echi.1 <- object1$chi

chi.2 <- object2$STATISTIC
dof.2  <- object2$dof
BIC.2  <- object2$BIC
echi.2 <- object2$chi

if(is.null(echi.1)) echi.1 <- NA
if(is.null(echi.2)) echi.2 <- NA
if( is.null(chi.1)) {stop("You do not seem to have chi square values for  ",name.1,"\nPerhaps you did not specify the sample size when you did the analysis?" )}
if( is.null(chi.2)) {stop("You do not seem to have chi square values for  ",name.2,"\nPerhaps you did not specify the sample size when you did the analysis?" )}
delta.df <- abs(dof.1 - dof.2)
delta.chi <-abs( chi.1 - chi.2)
delta.echi <- abs(echi.1 - echi.2)
test.chi <- delta.chi/delta.df
test.echi <- delta.echi/delta.df
delta.BIC <- (BIC.2 - BIC.1)




table <- data.frame(c(dof.1,dof.2),c(chi.1,chi.2),c(NA,delta.df),c(NA,delta.chi),c(NA, pchisq(delta.chi, delta.df, lower.tail=FALSE)),c(echi.1,echi.2),c(NA,delta.echi),c(NA,pchisq(delta.echi, delta.df, lower.tail=FALSE)),c(BIC.1,BIC.2),c(NA,delta.BIC))
names(table) <- c("Model Df", "ML Chisq", "Delta Df", "Delta Chisq", "Pr(> Delta Chisq)","Emp Chisq"," Delta Emp Chisq" ,"Pr(> Emp.Delta Chisq)","BIC","Delta BIC")
if(is.na(delta.echi)) {
table <- data.frame(c(dof.1,dof.2),c(chi.1,chi.2),c(NA,delta.df),c(NA,delta.chi),c(NA, pchisq(delta.chi, delta.df, lower.tail=FALSE)),c(BIC.1,BIC.2),c(NA,delta.BIC))
names(table) <- c("Model Df", "ML Chisq", "Delta Df", "Delta Chisq", "Pr(> Delta Chisq)","BIC","Delta BIC")
}
rownames(table) <- c(name.1, name.2)

structure(table,heading = c("ANOVA Test for Difference Between Models",""),
             class = c("anova", "data.frame"))		
			}
