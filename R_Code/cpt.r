cpt = function (Data,R,Crit) #[CP,L,Var]
{
# Uses t test to find first significant change point.
# Syntax is [CP,L,Var] = cpt(Data,R,Crit). Data is the vector of trial by
# trial measures or successive intervals. This test is appropriate if one
# is looking for a change in the expectation of a renewal event-generating
# process, where the interevent intervals are normally (rather than
# exponentially) distributed. Also when trial measures are normally
# distributed. R is the vector of putative change
# points. Crit is the decision criterion, the value the logit must exceed
# for the function to return a significant change point. CP is the
# first significant change point. L is the logit vector up to and including
# the row where the significance criterion (Crit) is exceeded. Var
# indicates the different possibilities for the variance estimates from the
# before and after data, to wit: Var==1 indicates that the estimates were
# both nonzero and not significantly different; Var==0 indicates that one
# of them was 0; Var==2 indicates that they were both nonzero and
# significantly different. The function uses the equal variance version of
# the t test in the Var=1 case (that is, it pools the variance
# estimates); it uses the unequal variance version of the t test
# in the second and third cases (Var = 0 0r 2). If there is no significant
# change point, CP is returned empty

  Cum = cumsum(Data) # The cumulative record
    
  L=matrix(0,nr=2,nc=1)
  r=3
  Var=matrix(c(0,0,0),nc=3)   # Initializing for while loop; r is the index
 
  while( abs(L[nrow(L),1])<Crit ) # loop that ends when critical L found or end of data reached
  {  
    Vb=var.ed(Data[1:R[r],1])
    dfb=length(Data[1:R[r],1])-1; # variance estimates
                                  # and associated degrees of freedom for the data up to and including
                                  # the putative change point
    
    Va= var.ed(Data[(R[r]+1):r,1])
    dfa=length(Data[(R[r]+1):r,1])-1 # ditto for data after putative change point   
        
    if (Vb+Va>0) # If #1: Is at least one estimate greater than zero?
    {            # test cannot be run when there is no variance on either side of putative CP, for
                 # example, in the sequence 0 0 0 3 3 3.
    
      if (min(Va,Vb)==0) # If #2: One estimate is zero, use unequal var t test
      {      
        if ((dfa*dfb) == 0)
        {
          r=r+1 
          next
        }                # Unequal variance version of t test requires that both df's be greater than 0, so when that
                         # condition is not satisfied, go on to next iteration
            
        pb = t.test(Data[1:R[r],1],Data[(R[r]+1):r,1],conf.level=.95,alternative = 'greater')$p.value
            # The .05 is a nominal significance level; it is there because
            # a third argument is required if there is to be a fourth. The
            # rest of the program uses only pb, not H.
            
        Var = rbind(Var,c(r, R[r], 0))
      } else if ( (pf(Va/Vb,dfa,dfb)>.95) || (pf(Va/Vb,dfa,dfb)<.05)) # The two non-zero variance estimates are significantly unequal
      {                                                        # (also use the unequal variance version of the t test). Note
                                                               # that for both variance estimates to be nonzero, both df's
                                                               # must also be nonzero
        pb = t.test(Data[1:R[r],1],Data[(R[r]+1):r,1],conf.level=.95,alternative = 'greater')$p.value
        Var = rbind(Var,c(r, R[r], 2)) # variance estimates both non-zero and significantly different
      } else # variances both nonzero and not significantly unequal Var=1; t test based on a single pooled variance estimate     
      {    
        pb = t.test(Data[1:R[r],1],Data[(R[r]+1):r,1],conf.level=.95,alternative = 'greater', var.equal = TRUE)$p.value 
      }                                     # equal variance version of if-elseif-else (If #2)

        L = rbind(L, log(pb/(1-pb), 10) ) #Latest logit
        
                
    }else
    {
      L = rbind(L,0)  # This else goes with If # 1 (are both variance estimates 0?)
    }                 # of If #1, which computes indvidual L values

    r=r+1; # Incrementing latest row for next iteration
    
    if (r>nrow(Data))  # terminate while loop if end of data reached
      break 
            
  } # of while loop; looping ceases either when end of data reached or when logit
    # exceeds decision criterion
        
  if ( abs( L[nrow(L),1] ) > Crit)
  {    
    CP=R[r-1,1];
  }  
  else
  {
    CP=NULL
  } 

  if ((var.ed(Data)>0) & (Va==0) & (Vb==0))
  {
    disp('Variance estimates up to and after change point both zero')
  }

return(CP)
}

var.ed = function(x)
{
  if(length(x)>1) {return(var(x))}
  else  return(0)
}