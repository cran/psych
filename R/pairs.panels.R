"pairs.panels" <-
function (x,smooth=TRUE,scale=FALSE,digits=2,pch=20,...) #combines a splom, histograms, and correlations
      { op <- par(no.readonly = TRUE) # the whole list of settable par's.
      par(pch=pch)
      if (smooth ){
         if (scale) {
             pairs(x,diag.panel=panel.hist,upper.panel=panel.cor.scale,lower.panel=panel.smooth,...)
                    }
                    else {
                    pairs(x,diag.panel=panel.hist,upper.panel=panel.cor,lower.panel=panel.smooth,...)
                          } 
                    }
                    else      #smooth is not true
             { if (scale) {pairs(x,diag.panel=panel.hist,upper.panel=panel.cor.scale,...)
               } else  {pairs(x,diag.panel=panel.hist,upper.panel=panel.cor,...) }
            } #end of else (smooth)
op <- par(op)
      }   #end of function

