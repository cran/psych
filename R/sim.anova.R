 "sim.anova" <-
function ( es1 = 0, es2 = 0, es3 = 0, es12 = 0, es13 = 0,
    es23 = 0, es123 = 0, es11=0,es22=0, es33=0,n = 2,n1 = 2, n2 = 2, n3 = 2, within=NULL,r=.8,factors=TRUE,center = TRUE,std=TRUE)
{
    contrasts <- function(n) {
        if (n%%2) {
            seq(-floor(n/2), floor(n)/2)
        }
        else {
            seq(-(n - 1), (n - 1), 2)
        }
    }

   if(n1 * n2 * n3) { n <- n * n1 * n2 * n3 }
   if(n1) { cont1 <- contrasts(n1)
           IV1 <- rep(cont1, n/n1)} else {IV1 <- IV1<- rnorm(n)}
   if(n2) {cont2 <- contrasts(n2)
           if (n1) {
           IV2 <- rep(outer(rep(1, n1), contrasts(n2)), n/(n2 * n1)) } else {
           IV2 <- rep(cont2,n/n2)}
           } else {IV2 <- rnorm(n)}
   if (n3) {cont3 <- contrasts(n3)
       if (n1) { if(n2) {
                    IV3 <- rep(outer(rep(1, n1 * n2), contrasts(n3)), n/(n1 *  n2 * n3))
                         } else {IV3 <- rep(outer(rep(1, n1 ), contrasts(n3)), n/(n1 *  n3))
                                 }
               } else {if(n2) {IV3 <- rep(outer(rep(1, n2 ), contrasts(n3)), n/(n2 *  n3)) } else {
                               IV3 <- rep(contrasts(n3),n/n3)}
        }
        } else {IV3=rnorm(n)}
        
     if(factors) {if(n1) {iv1<- factor(IV1)} else{iv1<- IV1}
               if(n2) {iv2<- factor(IV2)}  else{iv2<- IV2}
               if(n3) {iv3<- factor(IV3)} else{iv3<- IV3}
                }
    if(std) {IV1 <- IV1/sd(IV1)
            IV2 <-  IV2/sd(IV2)
            IV3 <-  IV3/sd(IV3)}
   
                         
       y <-  es1 * IV1 + es2 * IV2 + es3 * IV3 + es12 * IV1 * IV2 +
        es13 * IV1 * IV3 + es23 * IV2 * IV3 + es123 * IV1 * IV2 *
        IV3 + es11*IV1*IV1 + es22*IV2*IV2 + es33*IV3*IV3 + rnorm(n)
        
       if(!is.null(within)) {yw <- within
                             ny <- length(yw)
                         y<-t(r* t(matrix(rep(y,ny),ncol=ny))  + yw + sqrt(1-r^2) * rnorm(n))
                         } 
   

    if (!center) {
        IV1 <- IV1 - min(IV1) + 1
        IV2 <- IV2 - min(IV2) + 1
        IV3 <- IV3 - min(IV3) + 1
    }
    if(factors) {y.df <- data.frame(IV1=iv1,IV2=iv2,IV3=iv3,DV=y)} else {

    y.df <- data.frame(IV1, IV2, IV3,DV=y)}
}
