ystep = function(Cum)
{
  # YSTEP creates two-column array for plotting cumulative record when all y increments are 1
  # and x increments are real valued (usually, interevent intervals). Syntax
  # CR = ystep(Cum). Cum is a vector of cumulative interevent intervals.
  # The entry in the nth row of Cum is the interval from the onset of observation
  # to the occurrence of the nth event

  Cumd= matrix(c(0, rep(Cum,each = 2)))
  # makes doublet vector, which has two successive identical values for each single value of Cum
  # Puts t=0 at beginning of doublet time vector

  Nd=matrix(rep(c(1:length(Cum)) ,each = 2)) # Event-count vector
                                             # doublet vector for event counts

  Nd=rbind(0,0,matrix(Nd[-nrow(Nd)]) );
  # Count = 0 at t=0 and also at -t_1, the time immediately before the first event
  # Doublet count at end of record is dropped
 
  CR=cbind(Cumd, Nd); # Data required to make a step-plot cumulative record
  # with time (Cumd) of x axis and event count (Nd) on y axis

  return(CR)
}

