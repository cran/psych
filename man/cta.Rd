\name{cta}
\alias{cta}
\alias{cta.15}
\title{Simulate the C(ues) T(endency) A(ction) model of motivation}
\description{Dynamic motivational models such as the Dynamics of Action (Atkinson and Birch, 1970, Revelle, 1986) may be reparameterized as a simple pair of differential (matrix) equations (Revelle, 1986, 2008). This function simulates the dynamic aspects of the CTA.  The CTA model is discussed in detail in Revelle and Condon (2015).
}
\usage{
cta (n=3,t=5000, cues = NULL, act=NULL, inhibit=NULL,expect = NULL, consume = NULL, 
tendency = NULL,tstrength=NULL, type="both", fast=2,compare=FALSE,learn=TRUE,reward=NULL) 
cta.15(n = 3, t = 5000, cues = NULL,stim=NULL, act = NULL, inhibit = NULL, consume = NULL, 
   ten = NULL,  type = "both", fast = 2)


}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{number of actions to simuate }
  \item{t}{length of time to simulate}
  \item{cues}{a vector of cue strengths}
  \item{stim}{A vector of the environmental stimuli}
  \item{act}{matrix of associations between cues and action tendencies}
  \item{inhibit}{inhibition matrix }
  \item{consume}{ Consummation matrix }
  \item{ten}{Initial values of action tendencies }
  \item{type}{show actions, tendencies, both, or state diagrams }
  \item{fast}{display every fast time (skips  }
  \item{expect}{A matrix of expectations}
  \item{tendency}{starting values of tendencies}
  \item{tstrength}{a vector of starting value of tendencies}
  \item{compare}{Allows a two x two graph to compare two plots}
  \item{learn}{Allow the system to learn (self reinforce) over time}
  \item{reward}{The strength of the reward for doing an action}
}
\details{A very thorough discussion of the CTA model is available from Revelle (2008). An application of the model is discussed in Revelle and Condon (2015). 

\code{\link{cta.15}} is the version used to produce the figures and analysis in Revelle and Condon (2015).  \code{\link{cta}} is the most recent version and includes a learning function developed in collaboration with Luke Smillie at the University of Melbourne.

The dynamics of action (Atkinson and Birch, 1970)  was a model of how instigating forces elicited action tendencies which in turn elicited actions.  The basic concept was that action tendencies had inertia.  That is, a wish (action tendency) would persist until satisfied and would not change without an instigating force.  The consummatory strength of doing an action was thought  in turn to reduce the action tendency.  Forces could either be instigating or inhibitory (leading to "negaction").

 Perhaps the simplest example is the action tendency (T)  to eat a pizza.  The instigating forces (F) to eat the pizza include the smell and look of the pizza, and once eating it, the flavor and texture.  However, if eating the pizza, there is also a consummatory force (C) which was thought to reflect both the strength (gusto) of eating the pizza as well as some constant consummatory value of the activity (c).  If not eating the pizza, but in a pizza parlor, the smells and visual cues combine to increase the tendency to eat the pizza.  Once eating it, however, the consummatory effect is no longer zero, and the change in action tendency will be a function of both the instigating forces and the consummatory forces.  These will achieve a balance when instigating forces are equal to the consummatory forces.  The asymptotic strength of eating the pizza reflects this balance and does not require a ``set point" or ``comparator". 
 
 To avoid the problems of instigating and consummatory lags and the need for a decision mechanism, it is possible to reparameterize the original DOA model in terms of action tendencies and actions (Revelle, 1986).   Rather than specifying inertia for action tendencies and a choice rule of always expressing the dominant action tendency, it is useful to distinguish between action tendencies (t) and the actions (a) themselves and to  have actions as well as tendencies  having inertial properties. By separating tendencies from actions, and giving them both inertial properties, we avoid the necessity  of a  lag parameter, and by making the decision rule one of mutual inhibition, the process is perhaps easier to understand.  In an environment which affords cues for action (c), cues enhance action tendencies (t) which in turn strengthen actions (a).  This leads to two differential equations, one describing the growth and decay of action tendencies (t), the other of the actions themselves (a). 
 
 \deqn{d{t} = {Sc} - { Ca} }{dt = Sc - Ca}   and
 \deqn{d{a} = {Et} - {Ia}}{da = Et - Ia}.  
 (See Revelle and Condon (2015) for an extensive discussion of this model.)
 
 \code{\link{cta}} simulates this model, with the addition of a learning parameter such that  activities strengthen the connection between cues and tendencies.  The learning part of the cta model is still under development.  \code{\link{cta.15}} represents the state of the cta model as described in the Revelle and Condon (2015) article.
  
  



}
\value{
  graphical output unless type="none"
  
  \item{cues}{echo back the cue input}
  \item{inhibition}{echo back the inhibitory matrix}
  \item{time}{time spent in each activity}
  \item{frequency}{Frequency of each activity}
  \item{tendencies}{average tendency strengths}
  \item{actions}{average action strength}
  
}
\references{
Atkinson, John W. and Birch, David (1970) The dynamics of action. John Wiley, New York, N.Y.

Revelle, William (1986)  Motivation and efficiency of cognitive performance  in Brown, Donald R. and Veroff, Joe (ed). Frontiers of Motivational Psychology: Essays in honor of J. W. Atkinson. Springer.  (Available as a pdf at \url{https://personality-project.org/revelle/publications/dynamicsofmotivation.pdf}.)

Revelle, W. (2008)  Cues, Tendencies and Actions. The Dynamics of Action revisted. \url{https://personality-project.org/revelle/publications/cta.pdf} 

Revelle, W. and Condon, D. (2015) A model for personality at three levels.  Journal of Research in Personality, <doi:10.1016/j.jrp.2014.12.006>
}
\author{ William Revelle }


\examples{
#not run 
#cta()   #default values, running over time 
#cta(type="state") #default values, in a state space  of tendency 1 versus tendency 2
#these next are examples without graphic output
#not run
#two introverts
#c2i <- c(.95,1.05)
#cta(n=2,t=10000,cues=c2i,type="none")
#two extraverts
#c2e <- c(3.95,4.05)
#cta(n=2,t=10000,cues=c2e,type="none")
#three introverts
#c3i <-  c(.95,1,1.05)
#cta(3,t=10000,cues=c3i,type="none")
#three extraverts
#c3i <- c(3.95,4, 4.05)
#cta(3,10000,c3i,type="none")
#mixed
#c3 <- c(1,2.5,4)
#cta(3,10000,c3,type="none")
} 
\keyword{ models }

