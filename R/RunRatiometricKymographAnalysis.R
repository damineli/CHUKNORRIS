# Run ratiometric kymograph analysis scripts for a pair of kymographs
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
in1.path <- "~/Dropbox/CHUKNORRIS/data/ratio_kymo/"  # *** MUST SPECIFY! ***
in2.path <- "~/Dropbox/CHUKNORRIS/data/ratio_kymo/"  # *** MUST SPECIFY! ***
## Output folder - where do you want your data saved?
out.path <- "~/Dropbox/CHUKNORRIS/out_exs/ratio_kymo/"  # *** MUST SPECIFY! ***
# Files to analyze
## The files MUST be kymographs in '.txt', a matrix with pixel intensities where
## Time increases through rows (from top to bottom) and Object length increases
## through columns (from left to right)
################################################################################
# The simplified example below, with 1 being signal and 0 background, shows the
# correct orientation of the matrix for tip detection:
# 111111110000000000000000000000000000000000000000000000000000000000000000000000
# 111111111110000000000000000000000000000000000000000000000000000000000000000000
# 111111111111110000000000000000000000000000000000000000000000000000000000000000
# 111111111111111111100000000000000000000000000000000000000000000000000000000000
# 111111111111111111111111111111000000000000000000000000000000000000000000000000
# 111111111111111111111111111111111111111100000000000000000000000000000000000000
# 111111111111111111111111111111111111111111111111111100000000000000000000000000
################################################################################
## File name - what is the full name of the file 1 you want to analyze?
# File 1 MUST be the weakest channel, e.g. CFP
fl1.nm <- "Arab_Col-0_Kymo_020813-A_CFP.txt"  # *** MUST SPECIFY! ***
## Specify channel name in kymograph 1
var1.nm <- "CFP"  # can basically be anything between ""
# File 2 MUST be the strongest channel, e.g. YFP
fl2.nm <- "Arab_Col-0_Kymo_020813-A_YFP.txt"  # *** MUST SPECIFY! ***
## Specify variable name in kymograph 2
var2.nm <- "YFP"  # can basically be anything between ""
#-------------------------------------------------------------------------------
#--PARAMETERS (USER DEFINED)----------------------------------------------------
#-------------------------------------------------------------------------------
# Image parameters
image.par <- list(
                  # Time step between image acquisitions
                  time.step = 4 / 60,  #time step in minutes
                  # Unit of time in data 
                  time.unit = "min",   # currently either "min" or "s"
                  # Pixel size in length units
                  pixel.size = 0.2160014,  # here in micrometers
                  # Abbreviation of the unit used
                  length.unit = "um",  # micrometers
                  # Quantile of the background to subtract
                  bckgd.qntl = 0.99, # from 0 to 1, usually >> 0.75
                  # Number of pixels to use as a margin from the tip estimate
                  tip.mrgn = 0, # integer usually 0-5
                  # Number of pixels to measure Tip fluorescence counting from
                  # the tip estimate
                  tip.pixel = 10, # integer
                  # Number of pixels to measure Shank fluorescence counting from
                  # the tip estimate
                  shank.pixel = 90, # integer
                  # Number of pixels to extract fluorescence signal by averaging
                  # distributed symetrically around the tip and shank region
                  avg.width = 5, # integer, usually and odd number
                  # Number of pixels to cut from the tip for plotting purposes
                  tip.cut = 3 # integer < 5 to avoid problems with color scale
                  )

# Tip-finding parameters
tip.find.par <- list(
                     # Use smoothed fluoresce series to perform fit?
                     use.smooth = TRUE,  # TRUE or FALSE
                     # Number of points in the local fit to smooth kymograph
                     kymo.span.n = 7,  # integer ~3:21 (e.g. 7, 9)
                     # Degree of the polynomial used to smooth the kymograph
                     kymo.loess.dg = 2,  # integer 1, 2 (best 2)
                     # Number of points used in the linear model to find the tip
                     tip.find.n = 9,   # integer 5:11 (e.g. 7, 9, 11)
                     # Minumum number of points to consider a fluorescence 
                     # change (increase from right to left) chunk
                     fluo.chunk.n = 9,   # integer 7:11 (e.g. 7, 9, 11)
                     # Fraction of the fluorescence signal that must be reached
                     # withing the fluorescence chunk
                     fluo.thrsh.frac = 0.4,  #~0.1:0.75 (e.g. 0.4, 0.6)
                     # Number of points in the smoothing the tip location series
                     tip.span.n = 11,  # integer ~ 3:21 (e.g. 7, 9, 11)
                     # Degree of the polynomial used in tip location smoothing
                     tip.loess.dg = 2,   # integer 1,2
                     # Use smoothed tip series to align kymograph?
                     use.tip.smth = TRUE
                     )

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
  coarse.n = 6,  # integer usually around 5-20
  # Perform smoothing on an integrated version of the series?
  coarse.integrate = TRUE,  # TRUE or FALSE
  # Perform smoothing shifting all values to the positive range?
  coarse.positive = FALSE  # TRUE or FALSE
)

# Filter parameters
filter.par <- list(
  # Lowest period to analyze (minimum 2 * sampling rate)
  low.per = 8 / 60,  # in the same time units as data
  # Highest period to analyze (maximum 0.5 * total time span)
  high.per = 2,  # in the same time units as data
  # Vector of components to include or exclude in the discrete
  # wavelet tranform. Leave it NULL if you are clueless
  smth.vec = NULL # NULL or a vector with TRUE or FALSE as c(T, T, T, T, F, F)
)

# Wavelet parameters
wvlt.par <- list(
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

# Output parameters      
out.par <- list(
  # Do plotting?
  do.plot = TRUE,  # TRUE or FALSE
  # What label to use in the growth axis
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
#-------------------------------------------------------------------------------
#--OUTPUT DESCRIPTION-----------------------------------------------------------
#-------------------------------------------------------------------------------
# Four output tables will be saved: "OscillationAnalysisTbl" and "SummaryStats" 
# with all the time series generated in the module namely growth rate raw and 
# smoothed ('growth.raw' and 'growth.smth'), Tip and Shank fluorescence
# ('tip.fluo' and 'shank.fluo'), with explicit time series extracted presented
# in "RatioKymoTimeSeries" and the tip aligned ratiometric kymograph in 
# "TipAlignedKymo"

# Below DWT refers to Discrete Wavelet Transform 
# and CWT to Continuous Wavelet Transform

# OscillationAnalysisTbl" contain columns sequentially named: 
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

# "RatioKymoTimeSeries" contain columns sequentially named: 
# time - Just time in the original units
# tip.loc.raw - The raw estimates of tip location
# tip.loc.smth - Smoothed estimated of tip location   
# growth.raw - Growth rate calculated from raw tip estimates
# growth.smth - Growth rate calculated from smoothed tip estimates 
# tip.fluo - Tip fluorescence series 
# shank.fluo - Shank fluorescence series

# "TipAlignedKymo" is a matrix with ratio of fluorescence between both channels
# provided after backgroud subtraction.
# Time increases through rows (from top to bottom) and Object length increases
# through columns (from left to right) although the tip is aligned to the left
# with the base of the shank ending in the background to the right
################################################################################
# The simplified example below, with 1 being signal and 0 background, shows the
# orientation of the ratiometric kymograph matrix (same as for tip detection):
# 111111110000000000000000000000000000000000000000000000000000000000000000000000
# 111111111110000000000000000000000000000000000000000000000000000000000000000000
# 111111111111110000000000000000000000000000000000000000000000000000000000000000
# 111111111111111111100000000000000000000000000000000000000000000000000000000000
# 111111111111111111111111111111000000000000000000000000000000000000000000000000
# 111111111111111111111111111111111111111100000000000000000000000000000000000000
# 111111111111111111111111111111111111111111111111111100000000000000000000000000
################################################################################

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
requirements = c("biwavelet", "wavelets", "waveslim", "mixtools", "fields") #remove "waveslim"
LoadRequiredPackages(requirements = requirements)
#-------------------------------------------------------------------------------
#--READ FILE, RUN ANALYSIS FUNCTION AND SAVE OUTPUTS (AUTOMATIC)----------------
#-------------------------------------------------------------------------------
# Read kymograph 1
k1 <- as.matrix(read.table(paste(in1.path, fl1.nm, sep = ""), header = FALSE))
# Read kymograph 2
k2 <- as.matrix(read.table(paste(in2.path, fl2.nm, sep = ""), header = FALSE))

ratio.kymo.lst <- AnalyzeRatiometricKymograph(k1, k2, tip.find.par, image.par, 
                                              ts.par, filter.par, wvlt.par, 
                                              out.par, var1.nm, var2.nm,
                                              fl.nm = fl1.nm)
#-------------------------------------------------------------------------------