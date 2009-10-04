"alpha.scale" <-
function (x,y)   #find coefficient alpha given a scale and a data.frame of the items in the scale
	{
		n=length(y)          #number of variables
		Vi=sum(diag(var(y,na.rm=TRUE)))     #sum of item variance
		Vt=var(x,na.rm=TRUE)                #total test variance
		((Vt-Vi)/Vt)*(n/(n-1))}              #alpha

