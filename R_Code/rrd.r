rrd = function(Cum, R)
{

  # RRD computes the pseudologit vector for the discrete random rate case,
  # Syntax is  L=rrd(Cum,R)
  # Cum is the cumulative response count
  # R is the trial number of the putative change point
  # # L is the vector giving for each trial the pseudologit--log[P/(1-P+p)]--
  # for the probability of observing n_a or fewer post-CP responses if lambda_a/lambda_b = 1
  # where lambda_a and lambda_b are the rates of responding before and after
  # the putative change point; P is the probability of observing Cum-Cum(R) or fewer
  # responses after the putative change point; and p is the probability of
  # observing exactly Cum-Cum(R) responses after the putative change point.
  # Note that the numerator and denominator of the pseudologit are not complementary
  # probabilities, as they are in a true logit; they both include the probability of observing
  # exactly Cum-Cum(R) post-CP responses,
  # where Cum is the cumulative number of responsess at the moment of calculation.
  # The putative change point is the row (R=trial number) at which
  # the difference between the observed and expected number of responses is maximal
  # Cum/N is the slope of the cumulative record up to N (i.e., the average response rate
  # up to the trial of calculation
  # (Cum/N)*R is the expected number of responses up to trial R; Cum(R) is the observed number.
  # A negative pseudologit means that n_a is less than expected
 
  T=cbind(1:nrow(Cum)) # The trial count vector
 
  Ta=T-R # the trials since the putative change point
 
  p=Ta/T # probability of any one response occuring during the trials since the putative CP
 
  Na=Cum-Cum[R] # number of responses since putative CP

  Peqorl=pbinom(as.matrix(Na)[-1],as.matrix(Cum)[-1],p[-1]) 
  # Probability of observing Na or fewer total responses on the trials since putative CP
 
  Peqorm=1-Peqorl+dbinom(as.matrix(Na)[-1],as.matrix(Cum)[-1],p[-1]);#Probability of observing Na or more
  # total responses on the trials since putative CP. Note that this probability overlaps
  # the "complementary" prob; both include the probability of observing exactly Na events.
 
  L=cbind(c(0,log(Peqorl/Peqorm,10))); # Vector of the pseudologits
  return(L)
}