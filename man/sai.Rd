\name{sai}
\alias{sai}
\alias{tai}
\alias{sai.dictionary}
\docType{data}
\title{State Anxiety data from the PMC lab over multiple occasions.
}
\description{
State Anxiety was measured two-three times in 11 studies at the Personality-Motivation-Cognition laboratory.  Here are item responses for 11 studies (9 repeated twice, 2 repeated three times). In all studies, the first occasion was before a manipulation.  In some studies, caffeine, or movies or incentives were then given to some of the participants before the second and third STAI was given.  In addition, Trait measures are available and included in the tai data set (3032 subjects).  
}
\usage{data(sai)
data(tai)
data(sai.dictionary)
}
\format{
  A data frame with 3032 unique observations on the following 23 variables. 
  \describe{
    \item{\code{id}}{a numeric vector}
    \item{\code{study}}{a factor with levels \code{ages} \code{cart} \code{fast} \code{fiat} \code{film} \code{flat} \code{home} \code{pat} \code{rob} \code{salt} \code{shed}\code{shop} \code{xray}}
    \item{\code{time}}{1=First, 2 = Second, 3=third administration}
   \item{\code{TOD}}{TOD  (time of day 1= 8:50-9:30 am,2 = 1=3 pm, 3= 7:-8pm}
    \item{\code{drug}}{drug (placebo (0) vs. caffeine (1))}
     \item{\code{film}}{film (1=Frontline (concentration camp), 2 = Halloween
             3= National Geographic (control), 4- Parenthood (humor)}
    \item{\code{anxious}}{anxious}
    \item{\code{at.ease}}{at ease}
    \item{\code{calm}}{calm}
    \item{\code{comfortable}}{comfortable}
    \item{\code{confident}}{confident}
    \item{\code{content}}{content}
    \item{\code{high.strung}}{high.strung}
    \item{\code{jittery}}{jittery}
    \item{\code{joyful}}{joyful}
    \item{\code{nervous}}{nervous}
    \item{\code{pleasant}}{pleasant}
    \item{\code{rattled}}{over-excited and rattled}
    \item{\code{regretful}}{regretful}
    \item{\code{relaxed}}{relaxed}
    \item{\code{rested}}{rested}
    \item{\code{secure}}{secure}
    \item{\code{tense}}{tense}
    \item{\code{upset}}{upset}
    \item{\code{worried}}{worried}
    \item{\code{worrying}}{worrying}
  }
}
\details{The standard experimental study at the Personality, Motivation and Cognition (PMC) laboratory (Revelle and Anderson, 1997)  was to administer a number of personality trait and state measures (e.g. the \code{\link{epi}}, \code{\link{msq}}, \code{\link{msqR}} and \code{\link{sai}})  to participants before some experimental manipulation of arousal/effort/anxiety.  Following the manipulation (with a 30 minute delay if giving caffeine/placebo), some performance task was given, followed once again by measures of state arousal/effort/anxiety.  

Here are the item level data on the \code{\link{sai}} (state anxiety)  and the \code{\link{tai}} (trait anxiety). Scores on these scales may be found using the scoring keys.  The \code{\link{affect}} data set includes pre and post scores for two studies (flat and maps) which manipulated state by using four types of movies.

In addition to being useful for studies of motivational state, these studies provide examples of test-retest and  alternate form reliabilities.  Given that 10 items overlap with the \code{\link{msqR}} data, they also allow for a comparison of immediate duplication of items with 30 minute delays.

Studies CART, FAST, SHED,  RAFT,  and SHOP  were either control groups, or did not experimentally vary arousal/effort/anxiety.   

AGES, CITY, EMIT, RIM, SALT, and XRAY were caffeine manipulations between time 1 and 2
(RIM and VALE were repeated day 1 and day 2)

FIAT, FLAT, MAPS, MIXX, and THRU  were 1 day studies with film manipulation between time 1 and time 2.

SAM1 and SAM2 were the first and second day of a two day study. The STAI was given once per day.  MSQ not MSQR was given.

VALE and PAT were two day studies with the STAI given pre and post on both days

RIM was a two day study with the STAI and MSQ given once per day. 
 
Usually, time of day 1 = 8:50-9am am, and 2 = 7:30 pm, however, in rob, with paid subjects, the times were 0530 and 22:30.  

}
\source{Data collected at the Personality, Motivation, and Cognition Laboratory, Northwestern University, between 1991 and 1999.
}
\references{

Charles D. Spielberger and Richard L. Gorsuch and R. E. Lushene, (1970) Manual for the State-Trait Anxiety Inventory.

Revelle,  William  and Anderson, Kristen Joan  (1997) Personality, motivation and cognitive performance: Final report to the Army Research Institute on  contract MDA 903-93-K-0008

Rafaeli, Eshkol and Revelle, William (2006), A premature consensus: Are happiness and sadness truly opposite affects? Motivation and Emotion, 30, 1, 1-12.

Smillie, Luke D.  and Cooper, Andrew  and Wilt, Joshua  and Revelle, William (2012) Do Extraverts Get More Bang for the Buck? Refining the Affective-Reactivity Hypothesis of Extraversion. Journal of Personality and Social Psychology, 103 (2), 206-326.
	
}
\examples{
data(sai)
table(sai$study,sai$time)  #show the counts for repeated measures
#table(sai$study,sai$TOD) #and the studies by time of day
#table(sai$study,sai$drug) # and the studies with drug
#Here are the keys to score the sai total score, positive and negative items
sai.keys <- list(sai = c("tense","regretful" , "upset", "worrying", "anxious", "nervous" ,  
"jittery" , "high.strung", "worried" , "rattled","-calm", 
"-secure","-at.ease","-rested","-comfortable", "-confident" ,"-relaxed" , "-content" , 
"-joyful", "-pleasant"  ) ,
sai.p = c("calm","at.ease","rested","comfortable", "confident", "secure" ,"relaxed" ,     
       "content" , "joyful", "pleasant" ),  
sai.n = c( "tense" , "anxious", "nervous" , "jittery" , "rattled",     "high.strung",  
         "upset", "worrying","worried","regretful" )
) 

 #just get the control subjects
#control <- subset(sai,sai$study \%in\% c("Cart", "Fast", "SHED",  "RAFT", "SHOP")) 
#pre and post drug studies
#drug <- subset(sai,sai$study \%in\% c("AGES", "CITY","EMIT", "SALT", "VALE", "XRAY")) 
#pre and post film studies
#film <- subset(sai,sai$study \%in\% c("FIAT","FLAT", "MAPS", "MIXX") 

#this next set allows us to score those sai items that overlap with the msq item sets
msq.items <- c("anxious" ,  "at.ease" ,  "calm" ,     "confident", "content",   "jittery", 
 "nervous" ,  "relaxed" ,  "tense"  ,   "upset" ) #these overlap with the msq
 
sai.msq.keys <- list(pos =c( "at.ease" ,  "calm" , "confident", "content","relaxed"),
  neg = c("anxious", "jittery", "nervous" ,"tense"  ,   "upset"),
  anx = c("anxious", "jittery", "nervous" ,"tense", "upset","-at.ease" ,  "-calm" ,
  "-confident", "-content","-relaxed"))
sai.not.msq.keys <- list(pos=c(  "secure","rested","comfortable" ,"joyful" , "pleasant" ),    
    neg=c("regretful","worrying", "high.strung","worried", "rattled" ),
    anx = c("regretful","worrying", "high.strung","worried", "rattled",     "-secure",      
    "-rested", "-comfortable", "-joyful",  "-pleasant" )) 
sai.alternate.forms <- list( pos1 =c( "at.ease","calm","confident","content","relaxed"),
  neg1 = c("anxious", "jittery", "nervous" ,"tense"  ,   "upset"),
  anx1 = c("anxious", "jittery", "nervous" ,"tense", "upset","-at.ease" ,  "-calm" ,
       "-confident", "-content","-relaxed"),
  pos2=c(  "secure","rested","comfortable" ,"joyful" , "pleasant" ),    
  neg2=c("regretful","worrying", "high.strung","worried", "rattled" ),
  anx2 = c("regretful","worrying", "high.strung","worried", "rattled", "-secure",      
    "-rested", "-comfortable", "-joyful",  "-pleasant" )) 
  
#sai.repeated <- c("AGES","Cart","Fast","FIAT","FILM","FLAT","HOME","PAT","RIM","SALT",
#     "SAM","SHED","SHOP","VALE","XRAY")
#sai12 <- subset(sai,sai$study \%in\% sai.repeated)   #the subset with repeated measures
#Choose those studies with repeated measures by :
#sai.control <- subset(sai,sai$study %in% c("Cart", "Fast", "SHED", "SHOP"))
#sai.film <- subset(sai,sai$study %in% c("FIAT","FLAT") )  
#sai.drug <- subset(sai,sai$study %in% c("AGES",  "SALT", "VALE", "XRAY"))
#sai.day <- subset(sai,sai$study %in% c("SAM", "RIM"))
 
}
\keyword{datasets}
