#radar and spider developed July 10, 2011
#added the ability to rotate the entire figure 02/24/20 and to use keys
#Added label.pos to give better control of the figures
"radar" <- function(x,labels=NULL,keys=NULL,center=FALSE,connect=FALSE,scale=1,ncolors=31,fill=FALSE,add=FALSE,linetyp="solid", 
       main="Radar Plot",angle=0,absolute=FALSE,show=TRUE,digits=2,cut=.2,circles=TRUE,shape=FALSE, clockwise=FALSE,delta = NULL,label.pos=NULL,position=NULL,xlim=c(-1,1),ylim=c(-1, 1),
        ...) {
nvar <- length(x)
lty="solid"
 if(!is.null(keys) ) {key.ord <- selectFromKeys(keys)
 num.keys <- length(keys)
 for (i in 1:num.keys) {
 select <- sub("-", "", unlist(keys[i]))
 keys[[i]]<- select}
 key.ord <- unlist(keys)

} else {key.ord <-1:nvar
num.keys <- 1 
if (clockwise) key.ord <- nvar:1
               }

 if(is.null(labels)) {labels <- names(keys)
     show.keys<- TRUE} else {labels <- labels[key.ord]
     show.keys <- FALSE}


nvar <- length(key.ord)
x <- x[key.ord]  #this order the x vector according to the keys vector
if(is.null(labels)) labels <- paste("V",1:nvar,sep="")
if(is.null(delta)) {DELTA <- rep(1.1,nvar)} else {if(length(delta)==1) {DELTA =rep(delta,nvar) } else {DELTA=delta}}
if(is.null(label.pos)) {label.pos <- rep(1.05,nvar) } else {if(length(label.pos)==1) label.pos= rep(label.pos,nvar)}
if(is.null(position)) position <- rep(0,nvar)  
if (clockwise) { position <- position[key.ord]
                 label.pos <- label.pos[key.ord]
                 DELTA <- DELTA[key.ord]
                 
                }

SEGMENTS <- max(48,nvar)
if(ncolors < 2) {colors <- FALSE} else {colors <- TRUE}
 angles <- (0:SEGMENTS) * 2 * pi/SEGMENTS
 unit.circle <- cbind(cos(angles), sin(angles))

if(!add) {if(shape) {nx <- SEGMENTS/nvar
             Lx <- seq(1,SEGMENTS+1,nx)#%% SEGMENTS + 1
              plot(unit.circle[Lx,1],unit.circle[Lx,2],axes=FALSE,xlab="",ylab="",typ="l",xlim=xlim,ylim=ylim)

} else {plot(unit.circle,typ="l",asp=1,axes=FALSE,xlab="",ylab="",main=main,xlim=xlim,ylim=ylim)}
if(circles) {lines(unit.circle*.25,typ="l",lty="dotted",col="red")
lines(unit.circle*.5,typ="l",lty="dotted")
lines(unit.circle*.75,typ="l",lty="dotted",col="blue")} }

if(colors) { 
    gr <- colorRampPalette(c("red","white","blue")) #added June 20
    colramp  <- gr(ncolors)
      } else {
    colramp <- grey((ncolors:0)/ncolors)}

 for(c in 1:nvar) { 
  nx <- (((c-1 + angle*num.keys) %% nvar) +1 )* SEGMENTS/nvar +1
  nx <- nx %% SEGMENTS  # +1 
if(center) {x0 <- unit.circle[nx,1] * .5
            y0 <- unit.circle[nx,2] * .5
            } else {
           x0 <- 0
           y0 <- 0}

 if(!absolute) {scaler <- (x[c]*scale/2 + .5) } else {scaler <- abs(x[c]*scale    )}#stats can go from -1 to 1, scale from 0 to 1
 x1 <- unit.circle[nx,1]
 y1 <- unit.circle[nx,2] 
 Lx <- c(x0,x1)*scaler
 Ly <- c(y0,y1) *scaler
 if(c==1) {	Oldx <- unit.circle[(nvar-1)*SEGMENTS/nvar + 1,1]*(x[nvar]*scale/2 +.5)
 	        Oldy <- unit.circle[(nvar-1)*SEGMENTS/nvar + 1,2]*(x[nvar]*scale/2+.5)}
 
  if(colors) {   
   if(absolute) {
          {if (x[c] < 0 ) {col="red"
          lty="dotted"} else {col="blue"
          lty="solid"}
          }} else {if (scaler < .5) {col="red"} else {col="blue"}
          }
          

    
   lines(Lx,Ly,col=col,lty=lty,...) } else {
   lines(Lx,Ly,...)}
    if(show) {
        txt <- paste0(sprintf("%.2f",round(x[c],2)))
       if (abs(x[c]) > cut) {text(Lx[2]*DELTA[c],Ly[2]*DELTA[c],txt,...)}
       # if (abs(x[c]) > cut) {text(Lx[2]*sign(Lx[2])*DELTA[c],Ly[2]*sign(Ly[2])*DELTA[c],txt)}
      }     
   
 if(connect) {lines(c(Oldx,x1*scaler),c(Oldy,y1*scaler),lty=linetyp)}
 if(fill) {polygon(c(0,Oldx , x1*scaler,0),c(0,Oldy,y1*scaler,0),col=colramp[ceiling(scaler*ncolors)],...)} 
          
    Oldx <- x1*scaler
 	Oldy <- y1* scaler
 
 

if(!show.keys) { if(position[c] ==0) {pos=NULL} else { pos=position[c]}
          text(x1*label.pos[c],y1*label.pos[c],labels[c],pos=pos,...)} 
 if(show.keys) {for(c in 1:num.keys) {
    nx <- (((c-1 + angle) %% num.keys) +1 )* SEGMENTS/num.keys +1 
    x1 <- unit.circle[nx,1]
    y1 <- unit.circle[nx,2] 
    text(x1*label.pos[nx] ,y1*label.pos[nx],labels[c],...)
     }
 }
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
