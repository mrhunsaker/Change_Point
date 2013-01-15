slopes = function(CP)
{
# calculates x and y values for step plot of slopes
# Syntax S=slopes(CP)
# CP is a 2 column array, whose 1st col gives the row numbers of the
# successive significant change points and whose 2nd col gives
# the value of the cumulative count at that row number
# In other words, the input array gives the coordinates of the significant
# change points
# The 1st col of S contains the trial numbers at which there are
# significant changes in slopes.
# The 2nd col of S gives the slopes.
# In both cols, values are duplicated so as to make steps in the plotted
# function at the trials where there are significant change points

  if (nrow(CP)<2)
  {
    S=matrix(c(0,0,CP[1,1],0,CP[2,1]/CP[1,1],CP[2,1]/CP[1,1] ),nc=2)
  }else
  {
    Deltas = rbind(c(0,0), diff(CP)) #Delta trials (btwn Infl Pts) & delta Counts
    Slopes = as.matrix(Deltas[,2]/Deltas[,1]) #Successive slopes of the cum count record

    DblT=c(0, sort(rep(CP[,1],2))); # Trial (x) - axis doublet values

    DblT = as.matrix(DblT[-length(DblT)])  # Deletes 2nd half of last doublet trial

    Slopes1 =  rep(Slopes,each=2)
    S=cbind(DblT, Slopes1); # Returns array for plotting
    return(S)
  }
}
