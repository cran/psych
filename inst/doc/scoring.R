### R code from vignette source 'scoring.Rnw'

###################################################
### code chunk number 1: scoring.Rnw:164-166
###################################################
library(psych)
library(psychTools)


###################################################
### code chunk number 2: scoring.Rnw:206-208
###################################################
file.name <- "http://personality-project.org/r/psych/HowTo/scoring.tutorial/small.msq.txt"
my.data <- read.file(file.name)


###################################################
### code chunk number 3: scoring.Rnw:276-281
###################################################
fn <- "http://personality-project.org/r/psych/HowTo/scoring.tutorial/small.msq"
my.data <- read.file(fn, filetype="txt")  #because the suffix is not a standard one, we need to specify the file type
 dim(my.data)  #same as before
 headTail(my.data) #same as before
 describe(my.data)


###################################################
### code chunk number 4: scoring.Rnw:292-304
###################################################
 
my.keys <-  list(EA= c("active","alert","aroused", "-sleepy","-tired", "-drowsy"),
                            TA = c("anxious","jittery","nervous","-calm", "-relaxed", "-at.ease"),
                            EAp = c("active","alert","aroused"),
                            EAn = c("sleepy","tired", "drowsy"),
                            TAp = c("anxious","jittery","nervous"),
                            TAn = c("calm", "relaxed", "at.ease")
                            )
        
 another.keys.list <- list(EA=c(1:3,-4,-5,-6),TA=c(7:9,-10,-11,-12),
                    EAp =1:3,EAn=4:6,TAp =7:9,TAn=10:12)
   


###################################################
### code chunk number 5: scoring.Rnw:358-361
###################################################
my.scales <- scoreItems(my.keys,my.data)
my.scales   #show the output
my.scores <- my.scales$scores    #the actual scores are saved in the scores object


###################################################
### code chunk number 6: scoring.Rnw:416-417
###################################################
print(my.scales,short=FALSE)


###################################################
### code chunk number 7: scoring.Rnw:544-546
###################################################
my.scores <- my.scales$scores
headTail(round(my.scores,2) )


###################################################
### code chunk number 8: scoring.Rnw:554-556
###################################################
describe(my.scores)
pairs.panels(my.scores,pch='.')


###################################################
### code chunk number 9: pairs
###################################################
png('splom.png')
pairs.panels(my.scores,pch='.')
dev.off()


###################################################
### code chunk number 10: scoring.Rnw:587-593
###################################################
select <- colnames(my.data) 
#or
select <- selectFromKeys(my.keys)

small.msq <- msq[select]
describe(small.msq)  


###################################################
### code chunk number 11: scoring.Rnw:596-598
###################################################
msq.scales <- scoreItems(my.keys,small.msq)
msq.scales   #show the output


###################################################
### code chunk number 12: scoring.Rnw:659-661
###################################################
msq.scales.ov <- scoreOverlap(my.keys,small.msq)
msq.scales.ov   #show the output


###################################################
### code chunk number 13: scoring.Rnw:674-675
###################################################
sessionInfo()


