"VSS.sim" <-
function(ncases=1000,nvariables=16,nfactors=4,meanloading=.5,dichot=FALSE,cut=0)     #generates a simple structure factor matrix
                                                                    #with nfactors

{                                                                   
weight=sqrt(1-meanloading*meanloading)                            #loadings are path coefficients
theta=matrix(rnorm(ncases*nfactors),nrow=ncases,ncol=nvariables)  #generates nfactor independent columns, repeated nvariable/nfactor times)
error=matrix(rnorm(ncases*nvariables),nrow=ncases,ncol=nvariables) #errors for all variables
items=meanloading*theta+weight*error                               #observed score = factor score + error score
if(dichot) {items <- (items[,1:nvariables] >= cut) 
            items <- items + 0} 
return(items)
}


"VSS.simulate" <-
function(ncases=1000,nvariables=16,nfactors=4,meanloading=.5,dichot=FALSE,cut=0)     #generates a simple structure factor matrix
                                                                    #with nfactors

{                                                                   
weight=sqrt(1-meanloading*meanloading)                            #loadings are path coefficients
theta=matrix(rnorm(ncases*nfactors),nrow=ncases,ncol=nvariables)  #generates nfactor independent columns, repeated nvariable/nfactor times)
error=matrix(rnorm(ncases*nvariables),nrow=ncases,ncol=nvariables) #errors for all variables
items=meanloading*theta+weight*error                               #observed score = factor score + error score
if(dichot) {items <- (items[,1:nvariables] >= cut) 
            items <- items + 0} 
return(items)
}

