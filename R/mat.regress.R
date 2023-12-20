"mat.regress" <-
function(y,x,data,z=NULL,n.obs=NULL,use="pairwise",square=FALSE)  {
 #a function to extract subsets of variables (a and b) from a correlation matrix m or data set m
  #and find the multiple correlation beta weights + R2 of the a set predicting the b set
  #seriously rewritten, March 24, 2009 to make much simpler
  #minor additons, October, 20, 2009 to allow for print and summary function
  #major addition in April, 2011 to allow for set correlation
  
  message("mat.regress has been replaced by lmCor, please change your call") 
  lmCor(y,x,data,z=NULL,n.obs=NULL,use="pairwise",square=FALSE)} 

#modified July 12,2007 to allow for NA in the overall matrix
#modified July 9, 2008 to give statistical tests
#modified yet again August 15 , 2008 to convert covariances to correlations
#modified January 3, 2011 to work in the case of a single predictor 
#modified April 25, 2011 to add the set correlation (from Cohen)
