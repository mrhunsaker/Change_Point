flipud = function(diagonal)
 return(cbind(diagonal[nrow(diagonal):1,]))

cpd = function(Cum)
{
  #cpd finds putative change points in discrete-time cumulative records
  # It takes as input a vector of the cumulative response counts (or measures)
  # at ends of successive trials
  # Syntax is R=cpd(Cum)
  # The Nth row of Cum gives the count (or cumulative measure) up to and including the Nth trial
  # The putative change point corresponding to the Nth trial is the
  # preceding trial at which the deviation of the observed count or measure
  # from the expected count or measure is maximal. The expected count or measure
  # at end of an earlier trial, n, is the average count or measure per trial
  # over the range from n = 0 to n = N, times n.
  # The deviation from expectation is Cum(n) - this expectation. R is the value
  # of n at which this deviation is maximal
 
  N=cbind(1:nrow(Cum)) # Trial count vector
 
  Slopes=Cum/N  # Average count or measure per trial for trials 1 to N
 
  Dg=flipud(apply(flipud(diag(nrow(Cum))),1,cumsum))
  # Mask with ones on and above diagonal & zeros below
 
  Diagonal=Dg*t(matrix(t(Slopes),nrow(Slopes),nrow(Slopes)))
  # Creates an array in which successive cols have successive slopes of the cumulative record.
  # The slope for a given col fills all the cells on and above the main diagonal
 
  Preds=matrix(N,nrow(Cum),nrow(Cum))*Diagonal # Predicted (expected) cumulative values in a diagonal array
    
  Obs=matrix(t(Cum),nrow(Cum),nrow(Cum))*Dg; # Diagonal array of observed cumulations
 
  Devs=abs(Obs-Preds); # Diagonal array of deviations from expectations
 
  mx = apply(Devs,2, max ) # mx is a row vector listing the maximum in each col
  R = apply(Devs, 2, which.max) # R is a row vector specifying the row in which the max occurs
 
  R=cbind(R); # Converts R to col vector
  return(R)
}

