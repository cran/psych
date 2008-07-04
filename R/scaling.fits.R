"scaling.fits" <- 
function(model, data, test="logit", digits = 2,rowwise=TRUE) {
     model <- as.matrix(model)
     data <- as.matrix(data)
     if (test=="choice") {
     model <- as.vector(model)
     if (min(model) <= 0 ) model <- model -min(model)
     prob = model/(model %+% t(model)) 
     }  else {
     pdif <- model %+%-t(model)
      if (test=="logit") {
          prob <- 1/(1 + exp(-pdif))  }
     else {if (test=="normal") {
         prob <- pnorm(pdif)     }} }
      if (rowwise) {prob= 1- prob}
      error <- data - prob
  
     sum.error2 <- sum(error^2, na.rm = TRUE)
     sum.data2 <- sum(data^2, na.rm = TRUE)
     gof <- 1 - sum.error2/sum.data2
     fit <- list(GF = gof, original = sum.data2, resid = sum.error2,
     residual = round(error, 
          digits))
    return(fit)  }