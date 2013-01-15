rrc = function(Cum,R)
{
  #RRC computes pseudologit vector for continuous random rate case
  # Syntax is L=rrc(Cum,R)
  # Cum is the cumulative interevent interval vector
  # R is the vector specifying the row number (hence, the event number) of the putative CP
  # L is the vector giving for each event the pseudologit--log[P/(1-P+p)]--
  # for the probability of observing n_a or fewer post CP events if lambda_a/lambda_b = 1
  # where lambda_a and lambda_b are the rates of event occurrence before and after
  # the putative change point; P is the probability of observing N-R or fewer
  # events after the putative change point; and p is the probability of
  # observing exactly N-R events after the putative change point.
  # Note that the numerator and denominator of the pseudologit are not complementary
  # probabilities, as they are in a true logit; they both include the probability of observing
  # exactly N-R events, where N is the total number of events at the moment of calculation.
  # The putative change point is the row (R=event number) at which
  # the difference between (N/Cum(N))*(Cum(R) and Cum(R) is maximal
  # N/Cum(N) is the slope of the cumulative record up to the Nth event
  # Cum(R) is the time up to the Rth event
  # A negative pseudologit means that n_a is less than expected

  Ta=Cum-Cum[R] # the interval elapsed since the putative change point

  p=as.matrix(Ta/Cum) # probability of any one event falling after the putative change point

  N=cbind(1:nrow(Cum)); # event count vector

  Na=N-R; # vector giving number of events since putative change point

  Peqorl=pbinom(Na[-1],N[-1],p[-1])
  # Probability of observing Na or fewer events in the interval Ta

  Peqorm=1-Peqorl+dbinom(Na[-1],N[-1],p[-1]) #Probability of observing Na or more events
  # in the interval Ta. Note that this probability overlaps the "complementary" probability;
  # both include the probability of observing exactly Na events.

  L=cbind(c(0,log(Peqorm/Peqorl,10))) # Vector of the pseudologits
  return(L)
}
