trun =function (Cum,R,L,Crit)
{
  #TRUN truncates cumulative record at significant change points
  # Works for both cumulative responses vs trials and cumulative events vs time
  # that is, it works for both discrete and continuous cumulative records
  # Syntax is [Cumt Lt r]=trun(Cum,R,L,Crit)
  # All input arguments obligatory
  # Cum is cumulative vector (cumulative interevent intervals or cumulative responses)
  # R is putative change point vector
  # L is pseudologit vector, where logit is (approximately) the log of the odds that
  # there has been an CP;
  # Crit is decision criterion on logit;
  # Cumt is truncated cumulative record,
  # the record as it would be if observation began at time CP+ or after trial CP
  # Lt is the L vector truncated at Alert, which is the row at which an CP was detected
  # r is row at which it was truncated. All arguments are returned as empty
  # when there is no significant change point
 
  La=abs(L);

  Alert=min(which(La>Crit)) # Finds first row where decision criterion is exceeded
  # If there is no such row, then Alert will be empty

  if (Alert == Inf)
    Alert = NULL

  if (is.null(Alert)) #checks to see if Alert has any values
  { 
    Cumt=Cum; r=Alert; Lt=L
    return(list(Cumt,Lt,r)) ; # Returns the input Cum, input L & r = [] when there are no sig CPs
  }

  r=R[Alert] # The putative change point at the value of Cum that first yields significant logit.
  # This putative change point is always the number of an earlier event or trial
 
  I=matrix(c(Cum[1,],diff(as.matrix(Cum))) ) # Interevent interval vector OR vector of responses on successive trials
 
  Cumt=cumsum(I[(r+1):length(I)]); # Truncated cumulative record
 
  Lt=L[1:Alert] # Truncated logit vector
  list(Cumt,Lt,r)
}

