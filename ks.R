ks =function(Data,R,Crit) #[r1,r2,L,t]
{
  # KS computes the logit vector using the Kolmogorov-Smirnov test for
  # whether two distributions differ
  # The syntax is [r1,r2,L,t]=ks(Data,R,Crit)
  # Data is a column vector of successive measures (interevent intervals,
  # responses, poke durations, etc). NB, it is not the cumulative vector!
  # R is the vector of putative change points.
  # Crit is the decision criterion.
  # All arguments are obligatory.
  # L is the pseudologit vector for rows where the
  # approximation formula for the Kolmogorov-Smirnov p is valid.
  # Its final value is the first value to exceed the decision criterion
  # t is the col vector of rows for which L is defined
  # r1 is the change point row when the decision criterion is exceeded
  # r2 is the row at which the decision criterion is exceeded
  # If there are no testable rows or if the decision criterion is never
  # exceeded, the variables are returned empty
 
  N=as.matrix(c(1:nrow(Data))) # Number of rows in Data vector
 
  r1 = NULL
  r2 = NULL
  P = NULL
  L = NULL
  t = NULL # Initializing
 
  Na=N-R; # Col vector giving for each row in Data the number of rows
  #           after the putative change point. So R gives the number of
  #           rows before the change point and Na the number after
 
  Test=as.matrix(which((Na*R)/(Na+R)>=4)); # For the approximation to the KS
  #           probability to be valid, the product of the two n's divided by
  #           the sum must be greater than or equal to 4. Test is the col
  #           vector of rows that satisfy this constraint--the row numbers(!),
  #           not the entries themselves
 
  if (length(Test)==0) # There are no rows satisfying the constraint
   return(NULL) # Bail if there are no testable points
 
  for (T in 1:nrow(Test)) # Loop that steps through Data performing the KS test wherever
  {                       #   it is valid, until the resulting logit exceeds the decision
                          #   criterion
    P = rbind(P, ks.test(Data[1:R[Test[T]],1],Data[(R[Test[T]]+1):Test[T],1],alternative ='less')$p.value )
    L = rbind(L, abs( log( (1-P[T,1])/P[T,1],10 ) ) )
    t = cbind(t, Test[T,1]) # The row to which the latest value of L(T) "belongs"   

    if (L[T,1]>Crit) # Value of logit exceeds decision criterion
    {    
      r2=Test[T,1]; # the row (in Data) at which the criterion is exceeded
      r1=R[Test[T,1],1]; # the change point  when the criterion is exceeded,
      break # Break out of loop when decision criterion exceeded
    }# of if
    
  }# of for

  return(r1)
}

