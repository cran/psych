"panel.cor.scale" <-
function(x, y, digits=2, prefix="", cex.cor)
     {
         usr <- par("usr"); on.exit(par(usr))
         par(usr = c(0, 1, 0, 1))
         r = (cor(x, y,use="pairwise"))
         txt <- format(c(r, 0.123456789), digits=digits)[1]
         txt <- paste(prefix, txt, sep="")
         if(missing(cex.cor)) cex <- 0.8/strwidth(txt)
         cex1 <- cex  * abs(r)
         if(cex1 < .25) cex1 <- .25 #otherwise they just vanish
         text(0.5, 0.5, txt, cex = cex1)
     }

