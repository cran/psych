#radar and spider developed July 10, 2011
"radar" <- function(x,labels=NULL,center=FALSE,connect=FALSE,scale=1,ncolors=31,fill=FALSE,add=FALSE,linetyp="solid", main="Radar Plot",...) {
nvar <- length(x)
if(is.null(labels)) labels <- paste("V",1:nvar,sep="")

SEGMENTS <- 48
if(ncolors < 2) {colors <- FALSE} else {colors <- TRUE}
 angles <- (0:SEGMENTS) * 2 * pi/SEGMENTS
 unit.circle <- cbind(cos(angles), sin(angles))

if(!add) { plot(unit.circle,typ="l",asp=1,axes=FALSE,xlab="",ylab="",main=main)
 lines(unit.circle*.25,typ="l",lty="dotted",col="red")
lines(unit.circle*.5,typ="l",lty="dotted")
lines(unit.circle*.75,typ="l",lty="dotted",col="blue") }

if(colors) { 
    gr <- colorRampPalette(c("red","white","blue")) #added June 20
    colramp  <- gr(ncolors)
      } else {
    colramp <- grey((ncolors:0)/ncolors)}

 for(c in 1:nvar) {
  nx <-  (c-1)* SEGMENTS/nvar +1
if(center) {x0 <- unit.circle[nx,1] * .5
            y0 <- unit.circle[nx,2] * .5
            } else {
           x0 <- 0
           y0 <- 0}

 scaler <- (x[c]*scale/2 + .5)#stats can go from -1 to 1, scale from 0 to 1
 x1 <- unit.circle[nx,1]
 y1 <- unit.circle[nx,2] 
 Lx <- c(x0,x1)*scaler
 Ly <- c(y0,y1) *scaler
 if(c==1) {	Oldx <- unit.circle[(nvar-1)*SEGMENTS/nvar + 1,1]*(x[nvar]*scale/2 +.5)
 	        Oldy <- unit.circle[(nvar-1)*SEGMENTS/nvar + 1,2]*(x[nvar]*scale/2+.5)}
 
  if(colors) {   
   if (scaler < .5) {col="red"} else {col="blue"}
   lines(Lx,Ly,col=col,...) } else {
   lines(Lx,Ly,...)}

 if(connect) {lines(c(Oldx,x1*scaler),c(Oldy,y1*scaler),lty=linetyp)}
 if(fill) {polygon(c(0,Oldx , x1*scaler,0),c(0,Oldy,y1*scaler,0),col=colramp[ceiling(scaler*ncolors)],...)} 
          
    Oldx <- x1*scaler
 	Oldy <- y1* scaler
 	
text(x1*1.05,y1*1.05,labels[c])
 }
 }
 
"spider" <- function(y,x,data,labels=NULL,rescale = FALSE,center=FALSE,connect=TRUE,overlay=FALSE,scale=1,ncolors=31,fill=FALSE,main=NULL,...) {
if(is.null(labels))  labels <- colnames(data)[x]
if(rescale)  {
  data <- scale(data)/3 }   #rescales to -1 to 1
if(length(y)==1) {
    if(!is.null(main)) {main=main} else {main <- colnames(data)[y]}
    radar(data[y,x],labels=labels,center=center,connect=connect,scale=scale,ncolors=ncolors,fill=fill,main=main,...)  } else {
     nvar <- length(y)
for (i in 1:nvar) {
    if(!is.null(main)) {title=main[y[i]]} else {title <- colnames(data)[y[i]]}
  if(overlay) {
            if (i==1) {radar(data[y[i],x],labels=labels,center=center,connect=connect,scale=scale,ncolors=ncolors,fill=fill,main=title,...) } else {
                       radar(data[y[i],x],labels=labels,center=center,connect=connect,scale=scale,ncolors=ncolors,fill=fill,add=TRUE,linetyp=nvar %%6 + 2,main=title,...) }
              } else {
  
  radar(data[y[i],x],labels=labels,center=center,connect=connect,scale=scale,ncolors=ncolors,fill=fill,main=title,...) 
   }}
   }
}
