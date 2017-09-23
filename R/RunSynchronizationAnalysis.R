# Run synchronization analysis scripts for a pair of time series
################################################################################
# USER DEFINED SECTION #########################################################
################################################################################
#-------------------------------------------------------------------------------
#--INPUT/OUTPUT (USER DEFINED)--------------------------------------------------
#-------------------------------------------------------------------------------
# Directories 
## Path to CHUKNORRIS - where did you save the folder with the scripts?
chuk.path <- "~/Dropbox/"  # *** MUST SPECIFY! ***
## Input folders - where is your data?
in1.path <- "~/Dropbox/CHUKNORRIS/data/vp/"  # *** MUST SPECIFY! ***
in2.path <- "~/Dropbox/CHUKNORRIS/data/track/only_ts/"  # *** MUST SPECIFY! ***
## Output folder - where do you want your data saved?
out.path <- "~/Dropbox/CHUKNORRIS/out_exs/sync/"  # *** MUST SPECIFY! ***
# File to analyze
## The file MUST be a 2-column '.csv' file with Time in the first column and 
## variable values in the second column
## File name - what is the full name of the file 1 you want to analyze?
fl1.nm <- "Arab_Col-0_H-VP_190815-A.csv"  # *** MUST SPECIFY! ***
## Specify if the ".csv" file 1 is separated by commas "," or semicolon ";"
csv1.sep <- ";"  # options: "," or ";" *** MUST SPECIFY! *** 
## Specify variable name in file 1
var1.nm <- "H+_flux"  # can basically be anything between ""
## File name - what is the full name of the file 2 you want to analyze?
fl2.nm <- "TS_Object_1_Arab_Col-0_Track_190815-A.csv"  # *** MUST SPECIFY! ***
## Specify if the ".csv" file 2 is separated by commas "," or semicolon ";"
csv2.sep <- ","  # options: "," or ";" *** MUST SPECIFY! *** 
## Specify variable name in file 2
var2.nm <- "growth_rate"  # can basically be anything between ""
#-------------------------------------------------------------------------------
#--PARAMETERS (USER DEFINED)----------------------------------------------------
#-------------------------------------------------------------------------------
# There are 3 sets of parameters, being 2 for the individual analysis of 
# each series and 1 for the joint analysis. In this example the parameters
# for individual analysis are considered to be the same for simplicity, but one
# can easily defined them as separate lists

# Time series parameters for series 1
ts1.par <- list(
  # Unit of time in data 
  time.unit = "min",  # currently either "min" or "s"
  # Number of significant precision digits for time
  time.dgts = 2,  # an integer usually between 2-6
  # Use prefiltering module to interpolate/remove outliers/smooth
  do.prefilter = TRUE,  # do prefiltering? TRUE or FALSE
  # Number of points in the local fit to interpolate the series
  # your time series MUST have more than twice as many points
  interpol.n = 6,  # integer usually around 5-7
  # more for heavier smoothing 
  # Outlier removal technique to be used
  outlier.removal = "none",   # choices: 
  # "none", "tso", "mad", "trsh", "iqr"
  # Thershold for outliers removal if "thrsh" was chosen
  out.thrsh = NULL,  # has to be the absolute cut off value
  # Perform outlier removal on the integrated series? 
  outl.cum.ts = TRUE,  # TRUE or FALSE
  # Perform a coarser smoothing after interpolation?
  coarse.smth = TRUE,  # TRUE or FALSE
  # Number of points in the local fit to smooth the series
  coarse.n = 15,  # integer usually around 5-20
  # Perform smoothing on an integrated version of the series?
  coarse.integrate = TRUE,  # TRUE or FALSE
  # Perform smoothing shifting all values to the positive range?
  coarse.positive = FALSE  # TRUE or FALSE
)

# Time series parameters for series 2
# if desired to specify parameters independently, copy & paste the list above 
# replacing 'ts1.par' below and choose the desired values 
ts2.par <- ts1.par

# Filter parameters for series 1
filter1.par <- list(
  # Lowest period to analyze (minimum 2 * sampling rate)
  low.per = 16 / 60,  # in the same time units as data
  # Highest period to analyze (maximum 0.5 * total time span)
  high.per = 2,  # in the same time units as data
  # Vector of components to include or exclude in the discrete
  # wavelet tranform. Leave it NULL if you are clueless
  smth.vec = c(F, T, T, T, F, F),  # NULL or 
  # a vector with TRUE or FALSE
   # Use loess to remove trend
  loess.first = FALSE,
  # Smoothing span of the trend 
  loess.spn = 0.25,#0.4
   # Smoothing degree of the trend
  loess.degree = 1#0.2
)

# Filter parameters for series 2
# if desired to specify parameters independently, copy & paste the list above 
# replacing 'filter1.par' below and choose the desired values 
filter2.par <- filter1.par

# Wavelet parameters for pair analysis
wvlt.pair.par <- list(
  # Number of points in the local fit to interpolate both series
  # your time series MUST have more than twice as many points
  interpol.n = 6,  # integer usually around 5-7
  # more for heavier smoothing 
  # Significance level for detecting oscillatory components
  # corresponds to 1 - alpha
  p = 0.95, # standard = 1
  # How fine should be the period analysis?
  # this parameter controls the spacing between
  # wavelet scales
  dj = 1 / 100,  # the smaller the numBer, the finer resolution
  # Wavelet type as accepted by the 'biwavelet' package
  wvlt.type = "morlet",  # wavelet type
  # Pad the series with Zeros?
  do.pad = FALSE, # TRUE or FALSE, usually FALSE
  # What to filter out of the wavelet ridges, with the most
  # stringent being "coi", followed by "signif" or "none"
  filter.ridge = "coi", # "coi", "signif" or "none"
  # Maximum number of simulataneous periodic components to
  # write in the output table
  max.comp = 5,  # integer greater than 2, but hardly > 5 
  # Summarize wavelet ridges by fitting a mixture of normal
  # distributions?
  do.summary = FALSE, # TRUE or FALSE
  # Do wavelet coherence analysis for time series pairs
  do.wvlt.coher = FALSE,  # TRUE or FALSE
  # Number of randomization for wavelet coherence analysis
  nrands = 100  # integer usually between 100 and 1000
)

# Wavelet parameters for series 1 and 2
# if desired to specify parameters independently, copy & paste the list above 
# replacing 'wvlt.pair.par' below and choose the desired values 
wvlt1.par <- wvlt.pair.par
wvlt2.par <- wvlt.pair.par

# Output parameters for series 1              
out1.par <- list(
  # Do plotting?
  do.plot = TRUE,  # TRUE or FALSE
  # What label to use in the y axis
  ylab = expression(paste("H"^{"+"}," flux (pmol cm"^{-2},
                          " s"^{-1},")", sep = "")), # A string or
  # expression passable to 'ylab' in the base R plot function.
  # The simplest thing is to write something between quotation
  # marks (e.g. "m/s") without the need for special characters.
  # Save data tables outputs?     
  do.save.data = TRUE, # TRUE or FALSE
  # Where to save files?
  out.path = out.path # path between quotation marks or default
)

# Output parameters for series 2              
out2.par <- list(
  # Do plotting?
  do.plot = TRUE,  # TRUE or FALSE
  # What label to use in the y axis
  ylab = expression(paste("Growth rate (",mu, "m min"^{-1},
                          " s"^{-1},")", sep = "")), # A string or
  # expression passable to 'ylab' in the base R plot function.
  # The simplest thing is to write something between quotation
  # marks (e.g. "m/s") without the need for special characters.
  # Save data tables outputs?     
  do.save.data = TRUE, # TRUE or FALSE
  # Where to save files?
  out.path = out.path # path between quotation marks or default
)

# Group all parameters in a list to be passed to the function being executed
par.lst <- list(par1 = list(ts.par = ts1.par,
                            filter.par = filter1.par,
                            wvlt.par = wvlt1.par,
                            out.par = out1.par),
                par2 = list(ts.par = ts2.par,
                            filter.par = filter2.par,
                            wvlt.par = wvlt2.par,
                            out.par = out2.par))

# Output parameters for the joint analysis            
out.pair.par <- list(
  # Plot synchronization analysis outputs?  
  do.plot = TRUE,  # TRUE or FALSE
  # Save synchronization data tables outputs?        
  do.save.data = TRUE, # TRUE or FALSE
  # Where to save files?
  out.path = out.path # path between quotation marks or default
)
#-------------------------------------------------------------------------------
#--OUTPUT DESCRIPTION-----------------------------------------------------------
#-------------------------------------------------------------------------------
# Three output tables will be saved: "OscillationAnalysisTbl", "SyncTbl" and 
# "SummaryStats"
#
# Below DWT refers to Discrete Wavelet Transform 
# and CWT to Continuous Wavelet Transform
#
# "OscillationAnalysisTbl" contain columns sequentially named: 
# "fl.nm" - File name (user supplied)
# "var.nm" - Variable name (user supplied)
# "time" - Time after interpolation (if performed)
# "raw" - Raw data or interpolated/prefiltered data (if performed)
# "smth" - Data smoothed with DWT
# "osc" - Data smoothed and detrended with DWT (most relevant in the analysis)
# "trend" - Trend or baseline extracted with DWT
# "per" - Period of CWT ridge with maximum power
# "phs" - Phase of CWT ridge with maximum power
# "amp" - Amplitude of CWT ridge with maximum power
# "env" - Envelope of CWT ridge with maximum power
# "sig" - Significance of CWT ridge with maximum power
# "pwr" - Power of CWT ridge with maximum power
# "per.i" - Period of ith CWT ridge
# "phs.i" - Phase of ith CWT ridge
# "sig.i" - Significance of ith CWT ridge
# "pwr.i" - Power of ith CWT ridge
# "phs.tot" - Total phase summed of all i CWT ridges
# "amp.tot" - Total amplitude summed of all i CWT ridges
# "env.tot"- Total envelope summed of all i CWT ridges

# "SyncTbl" contain columns as "AnalysisTbl", differing in:
# fl.nm - Both files analyzed with "_VS_" in between
# var1.nm - Name of the series 1
# var2.nm - Name of the series 2
# technique - Joint wavelet analysis used 'cross.wvlt' (cross wavelet) or 
#                                         'wvlt.coher' (wavelet coherence)
# var1 - Values for the filtered (smoothed and detrended) and matched series 1
# var2 - Values for the filtered (smoothed and detrended) and matched series 2
# delay - Delay calculated for the wavelet ridge with maximum power
# delay.i - Delay calculated for the ith wavelet ridge

# "SummaryStats" contain rows for analysis performed for each series obtained in
# the filtering process with DWT, namely "raw", "smth", "osc" and "trend" (as
# described in the "AnalysisTbl" names above) with "osc" being the most relevant
# for the period analysis. The columns are sequentially named: 
# "fl.nm" - File name (user supplied)
# "var.nm" - Variable name (user supplied) with "_" and the DWT filtered series
# "min" - Minimum value
# "25%.quant" - The 25% quantile
# "median" - The median or 50% quantile
# "mean" - The mean
# "75%.quant" - The 75% quantile
# "max" - The maximum value
# "sd" - The standard deviation
# "mad" - The Median Absolute Deviation
# "iqr.1" - The Interquantile range between 25% and 75% quantile
# "1%.quant" - The 1% quantile
# "99%.quant" - The 99% quantile
# "iqr.2" - The Interquantile range between 1% and 99% quantile
# "period.pgram" - Fourier peak of the smoothed periodogram (most relevant)
# "period.ar" - Fourier peak of an AR model of the data

#-------------------------------------------------------------------------------
################################################################################
# AUTOMATIC SECTION ############################################################
################################################################################
#-------------------------------------------------------------------------------
#--LOAD ALL SCRIPTS (AUTOMATIC)-------------------------------------------------
#-------------------------------------------------------------------------------
src.path <- paste(chuk.path, "CHUKNORRIS/R/src/", sep = "")
all.scripts <- list.files(src.path, pattern = "[.]R$", recursive = TRUE)
invisible(lapply(paste(src.path, all.scripts, sep = ""), source))
#-------------------------------------------------------------------------------
#--CHECK, INSTALL AND LOAD REQUIRED PACKAGES (AUTOMATIC)------------------------
#-------------------------------------------------------------------------------
requirements = c("biwavelet", "wavelets", "mixtools") #remove "waveslim"
LoadRequiredPackages(requirements = requirements)
#-------------------------------------------------------------------------------
#--READ FILE, RUN ANALYSIS FUNCTION AND SAVE OUTPUTS (AUTOMATIC)----------------
#-------------------------------------------------------------------------------
# Read file1
dat1 <- read.csv(paste(in1.path, fl1.nm, sep = ""), sep = csv1.sep)
# Read file2
dat2 <- read.csv(paste(in2.path, fl2.nm, sep = ""), sep = csv2.sep)

joint.analysis.lst <- AnalyzeTimeSeriesSynchrony(dat1, dat2, par.lst, 
                                                 wvlt.pair.par, out.pair.par, 
                                                 var1.nm, var2.nm, 
                                                 fl1.nm, fl2.nm)
#-------------------------------------------------------------------------------