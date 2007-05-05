"paired.r" <-
function(xy,xz,yz,n) {
       diff <- xy-xz
       determin=1-xy*xy - xz*xz - yz*yz + 2*xy*xz*yz
       av=(xy+xz)/2
       cube= (1-yz)*(1-yz)*(1-yz)
       t2 = diff * sqrt((n-1)*(1+yz)/(((2*(n-1)/(n-3))*determin+av*av*cube)))
       return(t2)
        }

