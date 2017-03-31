# Run time series analysis scripts for a single time series
################################################################################
# USER DEFINED SECTION #########################################################
################################################################################
#-------------------------------------------------------------------------------
#--INPUT/OUTPUT (USER DEFINED)--------------------------------------------------
#-------------------------------------------------------------------------------
# Directories 
## Path to CHUKNORRIS - where did you save the folder with the scripts?
chuk.path <- "~/Dropbox/"  # *** MUST SPECIFY! ***
## Input folder - where is your data?
input.path <- "~/Dropbox/CHUKNORRIS/data/vp/"  # *** MUST SPECIFY! ***
## Output folder - where do you want your data saved?
out.path <- "~/Dropbox/CHUKNORRIS/out_exs/time_series/"  # *** MUST SPECIFY! ***
# File to analyze
## The file MUST be a 2-column '.csv' file with Time in the first column and 
## variable values in the second column
## File name - what is the full name of the file you want to analyze?
fl.nm <- "Arab_Col-0_H-VP_190815-A.csv"  # *** MUST SPECIFY! ***
## Specify if the ".csv" file is separated by commas "," or semicolon ";"
csv.sep <- ";"  # options: "," or ";" *** MUST SPECIFY! *** 
## Specify variable name
var.nm <- "H+_flux"  # can basically be anything between ""
#-------------------------------------------------------------------------------
#--PARAMETERS (USER DEFINED)----------------------------------------------------
#-------------------------------------------------------------------------------
# Time series parameters
ts.par <- list(
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

# Filter parameters
filter.par <- list(
                   # Lowest period to analyze (minimum 2 * sampling rate)
                   low.per = 16 / 60,  # in the same time units as data
                   # Highest period to analyze (maximum 0.5 * total time span)
                   high.per = 2,  # in the same time units as data
                   # Vector of components to include or exclude in the discrete
                   # wavelet tranform. Leave it NULL if you are clueless
                   smth.vec = c(F, T, T, T, F, F)  # NULL or 
                                                   # a vector with TRUE or FALSE
                   )

# Wavelet parameters
wvlt.par <- list(
                  # Significance level for detecting oscillatory components
                  # corresponds to 1 - alpha
                  p = 0.95, # standard = 1
                  # How fine should be the period analysis?
                  # this parameter controls the spacing between
                  # wavelet scales
                  dj = 1 / 100,  # the smaller the number, the finer resolution
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

# Output parameters               
out.par <- list(
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
#-------------------------------------------------------------------------------
#--OUTPUT DESCRIPTION-----------------------------------------------------------
#-------------------------------------------------------------------------------
# Two output tables will be saved: "OscillationAnalysisTbl" and "SummaryStats"
# Below DWT refers to Discrete Wavelet Transform 
# and CWT to Continuous Wavelet Transform

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
src.path <- paste(chuk.path, "/R/src/", sep = "")
all.scripts <- list.files(src.path, pattern = "[.]R$", recursive = TRUE)
invisible(lapply(paste(src.path, all.scripts, sep = ""), source))
#-------------------------------------------------------------------------------
#--CHECK, INSTALL AND LOAD REQUIRED PACKAGES (AUTOMATIC)------------------------
#-------------------------------------------------------------------------------
requirements = c("biwavelet", "wavelets", "mixtools") #remove "waveslim"
LoadRequiredPackages(requirements = requirements)
#-------------------------------------------------------------------------------
#--READ FILE AND RUN ANALYSIS FUNCTION (AUTOMATIC)------------------------------
#-------------------------------------------------------------------------------
# Read file
dat <- read.csv(paste(input.path, fl.nm, sep = ""), sep = csv.sep)

# Run analysis
analysis.lst <- AnalyzeTimeSeries(dat, ts.par, filter.par, wvlt.par, out.par, 
                                  var.nm = var.nm, fl.nm = fl.nm)
#-------------------------------------------------------------------------------
