#steps towards a Rasch modeling program for IRT 
#first, estimate the item difficulties 
"irt.item.diff.rasch" <-
function(items) {
 ncases <- nrow(items)
 item.mean <- colMeans(items,na.rm=TRUE)
 item.mean[item.mean<(1/ncases)] <- 1/ncases
 irt.item.diff.rasch <- log((1/item.mean)- 1) }
 