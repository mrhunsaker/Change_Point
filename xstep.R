xstep = function(Cum)
{
  #XSTEP creates two-column array for plotting cumulative record when x increment always 1
  # Syntax is CR=xstep(Cum)
  # Cum is a vector of cumulative response counts or measures over a sequence of trials
  # The entry in the nth row of Cum is the total count or measure as of end
  # of nth trial. This version of xstep corrects a defect in the original
  # version, which did not deal properly with cumulative records containing
  # negative values
 
  Nd=matrix(c(0,rep((1:length(Cum)),each=2))) # Duplicated trial-count vector in 1x2N array
                                              # This a column vector with each entry trial 
                                              # count a doublet entry. 0th trial is trial 
                                              #before observation begins; origin of x axis
 
  Cumd=matrix(rep(Cum,each=2))     # Duplicates and transposes Cum to make Colum vector 
                                   # with doublet entries for Cum
 
  Cumd=rbind(0,0,matrix(Cumd[-nrow(Cumd)]))  # There are 0 (observed) responses as of 0th trial
                            # and also at beginning of 1st trial. Cum(1) gives count as of end of 1st trial.
                            # Doublet of final cumulation is dropped
 
  CR=cbind(Nd, Cumd); # Data required to make a step-plot cumulative record
  # with trials (Nd) on x axis and cumulation (Cumd) on y axis

  return(CR)
}
