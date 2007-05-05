"VSS.simulate" <-
function(ncases,nvariables,nfactors,meanloading)     #generates a simple structure factor matrix
                                                                    #with nfactors

{                                                                   
weight=sqrt(1-meanloading*meanloading)                            #loadings are path coefficients
theta=matrix(rnorm(ncases*nfactors),nrow=ncases,ncol=nvariables)  #generates nfactor independent columns, repeated nvariable/nfactor times)
error=matrix(rnorm(ncases*nvariables),nrow=ncases,ncol=nvariables) #errors for all variables
items=meanloading*theta+weight*error                               #observed score = factor score + error score
return(items)
}

