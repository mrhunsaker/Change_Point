
#####################################
#Compiles all the programs together #
#####################################
source("cp_wrapper.r")
source("chi2logit.r")
source("cpc.r")
source("cpd.r")
source("cpt.r")
source("ks.r")
source("rdd.r")
source("rrc.r")
source("slopes.r")
source("trun.r")
source("xstep.r")
source("ystep.r")

#From now on, just run cp_wrapper(filename, isDiscrete, test, Crit)
#
# 
# cp_wrapper(filename,isDiscrete, test, Crit)
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




####################################################################################
#Reads the sample data from the online database / Not necessary to run the program #
####################################################################################
TestData1 = read.table("http://www.pnas.org/content/suppl/2004/08/31/0404965101.DC1/04965DataSet1.txt")
TestData2 = read.table("http://www.pnas.org/content/suppl/2004/08/31/0404965101.DC1/04965DataSet2.txt")
TestData3 = read.table("http://www.pnas.org/content/suppl/2004/08/31/0404965101.DC1/04965DataSet3.txt") #binary
TestData4 = read.table("http://www.pnas.org/content/suppl/2004/08/31/0404965101.DC1/04965DataSet4.txt")
TestData5 = read.table("http://www.pnas.org/content/suppl/2004/08/31/0404965101.DC1/04965DataSet5.txt") #binary
TestData6 = read.table("http://www.pnas.org/content/suppl/2004/08/31/0404965101.DC1/04965DataSet6.txt")

###############
#Sample Codes #
###############

cp_wrapper(TestData1, 0, 1, 2)
cp_wrapper(TestData2, 1, 3, 4)
cp_wrapper(TestData3, 1, 4, 2)
cp_wrapper(TestData4, 1, 2, 2)
cp_wrapper(TestData5, 1, 1, 3)
cp_wrapper(TestData6, 1, 2, 2)
cp_wrapper(TestData6, 1, 2, 1.3)
cp_wrapper(TestData6, 1, 3, 2)