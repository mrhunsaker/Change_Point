


cpc = function(Cum)
{
  #cpc finds putative change points in continuous-time cumulative records
  # It takes as input a vector of the cumulative interevent intervals
  # Syntax is  R=cpc(Cum)
  # The Nth row of Cum gives the interval from the onset of observation to
  # the Nth event. The putative change point corresponding to the Nth
  # event is the preceding event at which the deviation of the observed event count
  # from the expected event count is maximal. The expected event count
  # at any earlier event, n, is Cum(n), the interval up to the nth event,
  # divided by the average interevent interval over the range from n = 0 to n = N 
  # The deviation from expectation is n - this expectation. R is the value
  # of n at which this deviation is maximal
 
  N=cbind(1:nrow(Cum))
  # Event count vector
 
  Slopes=N/Cum  
  # Average slope up to given point in cumulative function
 
  Dg = flipud(apply(flipud(diag(nrow(Cum))),1,cumsum))

  Diagonal = Dg*t(matrix(t(Slopes), nrow(Slopes), nrow(Slopes)))
  # Creates an array in which successive cols have successive slopes of the cumulative record.
  # The slope for a given col fills all the cells on and above the main diagonal
 
  Preds=matrix(t(Cum),nrow(Cum),nrow(Cum))*Diagonal;
  # Creates diagonal array of the predicted numbers of events at each time in Cum
    
  Obs=matrix(N,nrow(Cum),nrow(Cum))*Dg; # Diagonal array with actual numbers of events
 
  Devs=abs(Obs-Preds); # Diagonal array of deviations from expectations
 
  mx = apply(Devs,2, max ) # mx is a row vector listing the maximum in each col
  R = apply(Devs, 2, which.max) # R is a row vector specifying the row in which the max occurs
 
  R=cbind(R) # Converts R to col vector
  return(R)
}

