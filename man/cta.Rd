\name{cta}
\alias{cta}
\title{Simulate the C(ues) T(endency) A(ction) model of motivation}
\description{Dynamic motivational models such as the Dynamics of Action (Atkinson and Birch, 1970, Revelle, 1986) may be reparameterized as a simple pair of differential (matrix) equations (Revelle, 1986, 2008). This function simulates the dynamic aspects of the CTA.  The CTA model is discussed in detail in Revelle and Condon (2015).
}
\usage{
cta(n = 3, t = 5000, cues = NULL, act = NULL, inhibit = NULL, consume = NULL, 
   ten = NULL,  type = "both", fast = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{number of actions to simuate }
  \item{t}{length of time to simulate}
  \item{cues}{a vector of cue strengths}
  \item{act}{matrix of associations between cues and action tendencies}
  \item{inhibit}{inhibition matrix }
  \item{consume}{ Consummation matrix }
  \item{ten}{Initial values of action tendencies }
  \item{type}{show actions, tendencies, both, or state diagrams }
  \item{fast}{display every fast time (skips  }

}
\details{A very thorough discussion of the CTA model is available from Revelle (2008). An application of the model is discussed in Revelle and Condon (2015). 
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

Revelle, W. (2008)  Cues, Tendencies and Actions. The Dynamics of Action revisted. \url{http://personality-project.org/revelle/publications/cta.pdf} 

Revelle, W. and Condon, D. (2015) A model for personality at three levels.  Journal of Research in Personality \url{http://personality-project.org/revelle/publications/jrp_14_final.pdf} 

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

