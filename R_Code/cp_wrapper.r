cp_wrapper = function(Data, isDiscrete, test, Crit)
{
# 
# CP=cp_wrapper(filename,isDiscrete, test, Crit)
# cp wrapper is used to combine the relevant functions for the changepoint
# algorithim into a single easy-to-use function. It returns an array giving
# the coordinates of the significant change points in the cumulative
# record. It graphs the cumulative record, with the change points
# indicated, and below that the slopes between the change points
#
# The first argument (the name of the file or variable containing the data) 
# is obligatory; it must be specified when the function is called.
# if the remaining three arguments are not supplied when the function is called,
# the user will be prompted for them. 
# If the first argument is a file to read data from, put single quotes around it.
# The file must be a tab-delimited ASCII file,
# with one datum per line, and with the .txt extension. If the first
# argument is a variable in a Matlab workspace, it must be a column vector. In
# either case, the data should be in uncumulated form (trial-by-trial
# measurements, or successive interevent intervals).
#
# Valid parameter values:
#   isDiscrete: 0 or 1
#     1 for discrete-trial data and 0 for successive real-valued intervals
#       (e.g., the intervals between successive events or the distances
#       between successive landmarks)
#   test: 1, 2, 3 or 4
#     1 for binomial (random rate) - must be integer-valued when used in
#        discrete-trials case
#     2 for Kolmogorov-Smirnov (real- or integer-valued data, although
#        technically only valid for real-valued variables)
#     3 for t test (either integer- or real-valued data)
#     4 for chisquare test (must be binary data - 0 or 1 in each row)
# 
#   Crit: A critical value use to test the significance of a given change point
#    Values should generally be in the range between 1.3, which corresponds to
#    a p value of .05 and 6, which corresponds to a p value of .000001



  #checks parameter
  if((length(match.call()) - 1) < 2) #Checks to see if only filename is input
  {
    isDiscrete = readline('Are these discrete-trial measurements? (1 if yes; 0 if no) \n') 
    test = readline('What test should be used to compare data before and after a putative change point? \n (Answer 1 for binomial, 2 for K-S, 3 for t test or 4 for chi square) \n')
    Crit = readline('Logit value (decision criterion between 1.3 and 6)? \n')
  }

  #consistency checks
  if(isDiscrete==0 & test==4)
  {
    cat('Cannot use chi square test when data are successive intervals \n')
    break
  }

  if (test==4 & sum(Data+!Data) != nrow(Data))
  {  
    cat('When chi square test is used, data must be binary, i.e. 0 or 1 \n')
    break
  }

  if (isDiscrete==1 & test==1 &  (sum(Data %% rep(1,nrow(Data))) != sum(rep(0,nrow(Data))))   )
  {
    cat('When the binomial test is used with discrete-trial data, the data must be integer valued')
    return
  }

  Cum=as.matrix(cumsum(Data)) #cumulative record

# Section for computing CP array when binomial (random rate) test is used
# or when chi square test is used


  if (test==1 | test==4)
  {  
    Cumt=Cum # Initializiing for while loop. The Cumt vector will
             # be truncated as change points are found
    
    if (test == 1)     # if binomial test is to be used
    {
        CritLength = 1 # When binomial test is used, there must be at least two data
    }else              # if chi square test is used
    {
      CritLength = 7 # A test of differences of frequency cannot be significant when
    }                # the total number of observations is less than 8
     
    CP=matrix(c(0, 0), nrow = 1)  # Initializing for while loop
    r=1    
    
    while (!is.null(r) & (nrow(Cumt)>CritLength))
    { 
      if (isDiscrete==0) # Data are continuous
      { 
          R=cpc(Cumt) # putative inflection points
          L=rrc(Cumt,R) # logit vector for continuous case
      }else # data are discrete
      {
          R=cpd(Cumt); # putative inflection points
          if (test==1) # if binomial test is to be used
          {    
            L=rrd(Cumt,R); # logit vector
          }else # if chisquare test is to be used
          {
            L=as.matrix(suppressWarnings(chi2logit(Cumt,R)))
          } # of computing logit vector in discrete case
      }  # of computing R & L for one pass

        Trunct = suppressWarnings(trun(Cumt,R,L,Crit)) #Cumt is the truncated cumulative
        Cumt = as.matrix(Trunct[[1]])   # record; Lt is the logit vector up to the point of
        Lt = Trunct[[2]]     # truncation (not used); r is the change point; r is empty if there is
        r = Trunct[[3]]      # no significant change point

         
        if (!is.null(r)) # if there is a change point, update change-point array
        {
          if (isDiscrete==0) #In the continuous case, the row count goes in the
          {                  # y-column of the output array (the event count); in all other cases, it goes
                             # in the x column. In the continuous case, the x column
                             # contains the successive event time
            CP =rbind(CP, c( Cum[CP[nrow(CP),2] + r,1 ] , CP[nrow(CP),2]+r))  # Add Cumt row for latest change point
                                                                              # to last change point to get Cum row of latest change point. 
                                                                              # Value of cumulative record at the change point
          }else # In the discrete case, the row data go in the first column of CP     
          {
            CP= rbind(CP, c( CP[nrow(CP),1]+r, Cum[CP[nrow(CP),1] +r ,1])) # Add Cumt row for latest change point
                                                                           # to last change point to get Cum row of latest change point
                                                                           # Value of cumulative record at the change point
          }
         } # of updating change-point array
         
    }# of while loop for finding successive change points when the binomial test is used
     
  }# of section that computes change-point array when binomial or chi square test is used




  if (test==2|test==3)
  {   
    NewData=Data; # Initializing for while loop. These vectors will
        # be truncated as change points are found
 
    CP=matrix(c(0, 0), nrow = 1)  # Initializing for while loop
    r=1    
   
    if (test==2)
    {
      CritLength=7; # K-S test is not valid when there
    }               # are fewer than 4 data in either of the two samples
    else
    {
      CritLength=2; # when t test is used there must be at least 3 data
    }
    
    while (!is.null(r) & (nrow(NewData)>CritLength) )
    {    
      NewCum=as.matrix(cumsum(NewData))
     
      if (isDiscrete==0) # Data are continuous
      {
        R=cpc(NewCum)
      }
      else # Data are discrete
      {
        R=cpd(NewCum)
      } # computing R
        
      if (test==2) # if K-S test is to be used
      {
        r=suppressWarnings(ks(NewData,R,Crit)) # r is the (significant) change point; if there is none, it's empty
      }
      else # if t test is to be used
      {
        r = cpt(NewData,R,Crit);
      }  # of computing new change point
 
      if (!is.null(r)) # if there is a new change point, update change-point array and truncate NewData             
      {
        CP = rbind(CP, c(CP[nrow(CP),1]+ r,Cum[CP[nrow(CP),1]+ r,1] ))  # Add Cumt row for latest change point
                                 # to last change point to get Cum row of latest change point
                                 # Value of cumulative record at the change point
        NewData=as.matrix(NewData[(r+1):nrow(NewData),1]) # Truncated data vector
             
      }# of updating change-point array & truncating
  
    }   # of while loop for computing CP array when K-S or t test are used
    
  }# of section that computes CP array when K-S or t test are used
 
# Adding final point to output array
  if (isDiscrete==1)
  {
    CP =rbind(CP, c(nrow(Cum), Cum[length(Cum),1]) ) # last row of CP array
  }                                                   # gives coordinates of final point in cumulative record
  else # in continuous case, row count goes in y column
  {
    CP = rbind(CP, c(Cum[nrow(Cum),1] ,nrow(Cum) ))
  } # of adding final point

#Computing Slopes

S = slopes(CP)

#Plotting
par(mfrow=c(2,1)) 

  if (isDiscrete == 0)
  {
    PlotArray = ystep(Cum)
  }else
  {
    PlotArray = xstep(Cum)
  }

  plot(PlotArray[,1], PlotArray[,2], t = 'l', ann = FALSE)
  points(CP[,1],CP[,2], cex  = 1.5)

  if  (isDiscrete == 0)
  {
    title(xlab = 'Time', ylab = 'Number of Events')
  }else
  {
    title(xlab = 'Trials', ylab = 'Cumulative response measure')
  }

  plot(S[,1],S[,2], t = 'l', ann = FALSE)

  if (isDiscrete == 0)
  {
    title(xlab = 'Time' , ylab = 'Average Rate = Events/(Unit Time)' )
  }else
  {
    title(xlab = 'Trials', ylab = 'Average Response per Trial')
  }
dimnames(CP)<-NULL
  return(CP)
}
