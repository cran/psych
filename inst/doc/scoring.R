## ----echo=TRUE----------------------------------------------------------------
library(psych)

## -----------------------------------------------------------------------------
file.name <-
  "https://personality-project.org/r/psych/HowTo/scoring.tutorial/small.msq.txt"
if(require(psychTools)) {my.data <- read.file(file.name)} else { my.data <- small.msq}

## -----------------------------------------------------------------------------
fn <- "https://personality-project.org/r/psych/HowTo/scoring.tutorial/small.msq"
if(require(psychTools)) {my.data <- read.file(fn, filetype="txt")} else {my.data <- small.msq}  #because the suffix is not a standard one, we need to specify the file type
 dim(my.data)  #same as before
 headTail(my.data) #same as before
 describe(my.data)

## ----echo=TRUE----------------------------------------------------------------
 
my.keys <-  list(EA= c("active","alert","aroused", "-sleepy","-tired", "-drowsy"),
            TA = c("anxious","jittery","nervous","-calm", "-relaxed", "-at.ease"),
            EAp = c("active","alert","aroused"),
            EAn = c("sleepy","tired", "drowsy"),
            TAp = c("anxious","jittery","nervous"),
            TAn = c("calm", "relaxed", "at.ease")
                            )
        
 another.keys.list <- list(EA=c(1:3,-4,-5,-6),TA=c(7:9,-10,-11,-12),
                    EAp =1:3,EAn=4:6,TAp =7:9,TAn=10:12)
   

## ----echo=TRUE----------------------------------------------------------------
my.scales <- scoreItems(my.keys,my.data)
my.scales   #show the output
my.scores <- my.scales$scores    #the actual scores are saved in the scores object

## ----echo=TRUE----------------------------------------------------------------
print(my.scales,short=FALSE)

## ----echo=TRUE----------------------------------------------------------------
my.scores <- my.scales$scores
headTail(round(my.scores,2) )

## ----echo=TRUE----------------------------------------------------------------
describe(my.scores)
pairs.panels(my.scores,pch='.')

## ----pairs, echo=FALSE--------------------------------------------------------
png('splom.png')
pairs.panels(my.scores,pch='.')
dev.off()

## ----echo=TRUE----------------------------------------------------------------
select <- colnames(my.data) 
#or

select <- selectFromKeys(my.keys)
if(require(psychTools)) {small.msq <- psychTools::msq[select]} else {small.msq <- my.data}
describe(small.msq) #note that if psychTools is not available this is just 200 cases


## ----echo=TRUE----------------------------------------------------------------
msq.scales <- scoreItems(my.keys,small.msq)
msq.scales   #show the output

## ----echo=TRUE----------------------------------------------------------------
msq.scales.ov <- scoreOverlap(my.keys,small.msq)
msq.scales.ov   #show the output

## ----echo=TRUE----------------------------------------------------------------
sessionInfo()

