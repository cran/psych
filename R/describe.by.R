"describe.by" <-
function (x,group,...) {               #data are x, grouping variable is group
answer <- by(x,group,describe,...)
return(answer)}