"circ.sim.plot" <- 
  function(x.df) { 
  with(x.df,{
  symb <- c(21,22,20)
  colors <- c("black","blue", "red")
 op <- par(mfrow=c(2,2))
  plot(c.gap,c.RT,xlim=c(0,.5),ylim=c(0,1),pch=symb[1],col=colors[1],xlab="Gap Test",ylab="Rotation Test",main="Gap x Rotation")
  points(s.gap,s.RT,pch=symb[2],col=colors[2])
    points(e.gap,e.RT,pch=symb[3],col=colors[3])
  
 plot(c.gap,c.fisher,xlim=c(0,.5),ylim=c(0,.50),pch=symb[1],col=colors[1],xlab="Gap Test",ylab="Fisher Test",main="Gap x Fisher")
  points(s.gap,s.fisher,pch=symb[2],col=colors[2])
   points(e.gap,e.fisher,pch=symb[3],col=colors[3])
  
   plot(c.fisher,c.RT,xlim=c(0,.5),ylim=c(0,1),pch=symb[1],col=colors[1],xlab="Fisher Test",ylab="Rotation Test",main="Fisher x Rotation")
  points(s.gap,s.RT,pch=symb[2],col=colors[2])
   points(e.gap,e.RT,pch=symb[3],col=colors[3])
   
   boxplot(x.df,main=" Box Plot of all tests")
  title(main = "Circumplex Tests for Circumplex, Ellipsoid, and Simple Structure",outer=TRUE,line=-1)
})  #end of with
op <- par(mfrow=c(1,1))
  }