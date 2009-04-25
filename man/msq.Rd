\name{msq}
\Rdversion{1.1}
\alias{msq}
\docType{data}
\title{75 mood items from the Motivational State Questionnaire for 3896 participants}
\description{Emotions may be described either as discrete emotions or in dimensional terms.  The Motivational State Questionnaire (MSQ) was developed to study emotions in laboratory and field settings.  The data can be well described in terms of a two dimensional solution of energy vs tiredness and tension versus calmness.  Additional items include what time of day the data were collected and a few personality questionnaire scores.
}
\usage{data(msq)}
\format{
  A data frame with 3896 observations on the following 92 variables.
  \describe{
    \item{\code{MSQ_Time}}{Time of day the data were collected}
    \item{\code{active}}{a numeric vector}
    \item{\code{afraid}}{a numeric vector}
    \item{\code{alert}}{a numeric vector}
    \item{\code{alone}}{a numeric vector}
    \item{\code{angry}}{a numeric vector}
    \item{\code{aroused}}{a numeric vector}
    \item{\code{ashamed}}{a numeric vector}
    \item{\code{astonished}}{a numeric vector}
    \item{\code{at-ease}}{a numeric vector}
    \item{\code{at-rest}}{a numeric vector}
    \item{\code{attentive}}{a numeric vector}
    \item{\code{blue}}{a numeric vector}
    \item{\code{bored}}{a numeric vector}
    \item{\code{calm}}{a numeric vector}
    \item{\code{clutched-up}}{a numeric vector}
    \item{\code{confident}}{a numeric vector}
    \item{\code{content}}{a numeric vector}
    \item{\code{delighted}}{a numeric vector}
    \item{\code{depressed}}{a numeric vector}
    \item{\code{determined}}{a numeric vector}
    \item{\code{distressed}}{a numeric vector}
    \item{\code{drowsy}}{a numeric vector}
    \item{\code{dull}}{a numeric vector}
    \item{\code{elated}}{a numeric vector}
    \item{\code{energetic}}{a numeric vector}
    \item{\code{enthusiastic}}{a numeric vector}
    \item{\code{excited}}{a numeric vector}
    \item{\code{fearful}}{a numeric vector}
    \item{\code{frustrated}}{a numeric vector}
    \item{\code{full-of-pep}}{a numeric vector}
    \item{\code{gloomy}}{a numeric vector}
    \item{\code{grouchy}}{a numeric vector}
    \item{\code{guilty}}{a numeric vector}
    \item{\code{happy}}{a numeric vector}
    \item{\code{hostile}}{a numeric vector}
    \item{\code{inspired}}{a numeric vector}
    \item{\code{intense}}{a numeric vector}
    \item{\code{interested}}{a numeric vector}
    \item{\code{irritable}}{a numeric vector}
    \item{\code{jittery}}{a numeric vector}
    \item{\code{kindly}}{a numeric vector}
    \item{\code{lively}}{a numeric vector}
    \item{\code{lonely}}{a numeric vector}
    \item{\code{nervous}}{a numeric vector}
    \item{\code{placid}}{a numeric vector}
    \item{\code{pleased}}{a numeric vector}
    \item{\code{proud}}{a numeric vector}
    \item{\code{quiescent}}{a numeric vector}
    \item{\code{quiet}}{a numeric vector}
    \item{\code{relaxed}}{a numeric vector}
    \item{\code{sad}}{a numeric vector}
    \item{\code{satisfied}}{a numeric vector}
    \item{\code{scared}}{a numeric vector}
    \item{\code{scornful}}{a numeric vector}
    \item{\code{serene}}{a numeric vector}
    \item{\code{sleepy}}{a numeric vector}
    \item{\code{sluggish}}{a numeric vector}
    \item{\code{sociable}}{a numeric vector}
    \item{\code{sorry}}{a numeric vector}
    \item{\code{still}}{a numeric vector}
    \item{\code{strong}}{a numeric vector}
    \item{\code{surprised}}{a numeric vector}
    \item{\code{tense}}{a numeric vector}
    \item{\code{tired}}{a numeric vector}
    \item{\code{unhappy}}{a numeric vector}
    \item{\code{upset}}{a numeric vector}
    \item{\code{vigorous}}{a numeric vector}
    \item{\code{wakeful}}{a numeric vector}
    \item{\code{warmhearted}}{a numeric vector}
    \item{\code{wide-awake}}{a numeric vector}
    \item{\code{anxious}}{a numeric vector}
    \item{\code{idle}}{a numeric vector}
    \item{\code{cheerful}}{a numeric vector}
    \item{\code{inactive}}{a numeric vector}
    \item{\code{tranquil}}{a numeric vector}
    \item{\code{EA}}{Thayer's Energetic Arousal Scale}
    \item{\code{TA}}{Thayer's Tense Arousal Scale}
    \item{\code{PA}}{Positive Affect scale}
    \item{\code{NegAff}}{Negative Affect scale}
    \item{\code{Extraversion}}{Extraversion from the Eysenck Personality Inventory}
    \item{\code{Neuroticism}}{Neuroticism from the Eysenck Personality Inventory}
    \item{\code{Lie}}{Lie from the EPI}
    \item{\code{Sociability}}{The sociability subset of the Extraversion Scale}
    \item{\code{Impulsivity}}{The impulsivity subset of the Extraversions Scale}
    \item{\code{MSQ_Round}}{Rounded time of day}
    \item{\code{scale}}{a factor with levels \code{msq} \code{r} original or revised msq}
    \item{\code{ID}}{subject ID}
    \item{\code{exper}}{Which study were the data collected: a factor with levels \code{AGES} \code{BING} \code{BORN} \code{CART} \code{CITY} \code{COPE} \code{EMIT} \code{FAST} \code{Fern} \code{FILM} \code{FLAT} \code{Gray} \code{imps} \code{item} \code{knob} \code{MAPS} \code{mite} \code{pat-1} \code{pat-2} \code{PATS} \code{post} \code{RAFT} \code{Rim.1} \code{Rim.2} \code{rob-1} \code{rob-2} \code{ROG1} \code{ROG2} \code{SALT} \code{sam-1} \code{sam-2} \code{SAVE/PATS} \code{sett} \code{swam} \code{swam-2} \code{TIME} \code{VALE-1} \code{VALE-2} \code{VIEW}}
    \item{\code{condition}}{What was the experimental condition after the msq was given}
    \item{\code{TOD}}{a numeric vector}
    \item{\code{TOD24}}{a numeric vector}
  }
}
\details{The Motivational States Questionnaire (MSQ) is composed of 72 items, which represent the full affective range (Revelle & Anderson, 1996). The MSQ consists of 20 items taken from the Activation-Deactivation Adjective Check List (Thayer, 1986), 18 from the Positive and Negative Affect Schedule (PANAS, Watson, Clark, & Tellegen, 1988) along with the items used by Larsen and Diener (1992). The response format was a four-point scale that corresponds to Russell and Carroll's (1999) "ambiguous--likely-unipolar format" and that asks the respondents to indicate their current standing (``at this moment") with the following rating scale:\cr
0----------------1----------------2----------------3
\cr
Not at all		A little		Moderately     	Very much \cr
The original version of the MSQ included 70 items. Intermediate analyses (done with 1840 subjects) demonstrated a concentration of items in some sections of the two dimensional space, and a paucity of items in others. To begin correcting this, 3 items from redundantly measured sections (alone, kindly, scornful) were removed, and 5 new ones (anxious, cheerful, idle, inactive, and tranquil) were added. 

Procedure. The data were collected over nine years, as part of a series of studies examining the effects of personality and situational factors on motivational state and subsequent cognitive performance. In each of 38 studies, prior to any manipulation of motivational state, participants signed a consent form and filled out the MSQ. (The procedures of the individual studies are irrelevant to this data set and could not affect the responses to the MSQ, since this instrument was completed before any further instructions or tasks).

}
\source{Data collecte at the Personality, Motivation, and Cognition Laboratory, Northwestern University. 
}
\references{

William Revelle and Kristen Joan Anderson (1997) Personality, motivation and cognitive performance: Final report to the Army Research Institute on  contract MDA 903-93-K-0008

Rafaeli, Eshkol and Revelle, William (2006), A premature consensus: Are happiness and sadness truly opposite affects? Motivation and Emotion, 30, 1, 1-12.
}	

\examples{
data(msq)
describe(msq)
}
\keyword{datasets}
