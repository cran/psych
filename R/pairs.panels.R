"pairs.panels" <-
function (x,y,smooth=TRUE,scale=FALSE,digits=2,...) #combines a splom, histograms, and correlations
      {if (smooth ){
         if (scale) {
             pairs(x,diag.panel=panel.hist,upper.panel=panel.cor.scale,lower.panel=panel.smooth,...)
                    }
                    else {pairs(x,diag.panel=panel.hist,upper.panel=panel.cor,lower.panel=panel.smooth,...)
                    } #else  {pairs(x,diag.panel=panel.hist,upper.panel=panel.cor,lower.panel=panel.smooth,...)
                    }
                   
                    else      #smooth is not true
             { if (scale) {pairs(x,diag.panel=panel.hist,upper.panel=panel.cor.scale,...)
               } else  {pairs(x,diag.panel=panel.hist,upper.panel=panel.cor,...) }
            } #end of else (smooth)
         
      }   #end of function

