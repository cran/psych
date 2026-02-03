#find the Effective dimensionality of data set using the DelGuidice n1 statistic
#Developed 12/30/25
ED <- function(r,cor="cor", use="pairwise", n.obs=NA,weight=NULL,correct=.0) {


 if(!isCorrelation(r)) {n.obs <- nrow(r)
 original.data <- r   #save these for scores 
  switch(cor, 
       cor = {if(!is.null(weight))  {r <- cor.wt(r,w=weight)$r} else  {
                                     r <- cor(r,use=use)}
                                     },
       cov = {if(!is.null(weight))  {r <- cor.wt(r,w=weight,cor=FALSE)$r} else  {
                                     r <- cov(r,use=use)}
                                       covar <- TRUE},  #fixed 10/30/25
       wtd = { r <- cor.wt(r,w=weight)$r},
       spearman = {r <- cor(r,use=use,method="spearman")},
       kendall = {r <- cor(r,use=use,method="kendall")},
       tet = {r <- tetrachoric(r,correct=correct,weight=weight)$rho},
       poly = {r <- polychoric(r,correct=correct,weight=weight)$rho},
       tetrachoric = {r <- tetrachoric(r,correct=correct,weight=weight)$rho},
       polychoric = {r <- polychoric(r,correct=correct,weight=weight)$rho},
       mixed = {r <- mixedCor(r,use=use,correct=correct)$rho},
       Yuleb = {r <- YuleCor(r,,bonett=TRUE)$rho},
       YuleQ = {r <- YuleCor(r,1)$rho},
       YuleY = {r <- YuleCor(r,.5)$rho } 
       )
       
     R <- S <- r
     
         } else {R <- S <- r
     if(is.na(n.obs)){ n.obs <- 100
         cat("\n Number of observations not specified. Arbitrarily set to 100.\n")}
     original.data  <- r}
     
    ev <- eigen(r)
    ev.values <- ev$values
    tot <- sum(ev.values)
    
    
    nfactors<- 8
     fa.valueso  <- fa(r,nfactors=nfactors, warnings=FALSE)$Vaccounted[1,]
     fa.valuesn  <- fa(r,nfactors=nfactors,rotate="none",warnings=FALSE)$Vaccounted[1,]
     fa.valuesv <-  fa(r,nfactors=nfactors,rotate="simplimax", warnings=FALSE)$Vaccounted[1,]
     tot.values <- sum(fa.valueso)
     tot.valuesv <- sum(fa.valuesv)
     tot.valuesn <- sum(fa.valuesn)
    p <- ev.values/tot    #p as eigen value converted to a probability
    pob <-fa.valueso/tot.values
    pno <- fa.valuesn/tot.valuesn
    pva  <- fa.valuesv/tot.valuesv
    Hno <- -sum(pob*log(pob))
     Hnn <- -sum(pno*log(pno))
     Hnv <- -sum(pva*log(pva))
     
     pc <- pca(r,nfactors=nfactors)$Vaccounted[1,]
     pcv <- sum(pc)
     pc <- pc/pcv
     Hpc <- -sum(pc * log(pc))
    H <- -sum(p * log(p))   #shannon coefficient 
    n1 <- exp(H)
    redundancy <- 1-n1/tot
    result<- list(n1=n1,redundancy=redundancy,totob=tot.values,totvar=tot.valuesv,totnon=tot.valuesn, Hno=Hno,Hnv=Hnv,Hnn=Hnn,nobl=exp(Hno),nsiml=exp(Hnv),nnon=exp(Hnn),npc = exp(Hpc))
    return(result)
 }