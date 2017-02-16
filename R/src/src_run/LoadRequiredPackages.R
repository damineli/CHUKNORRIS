#-------------------------------------------------------------------------------
LoadRequiredPackages <- function(requirements = c("biwavelet", "wavelets", 
                                                  "waveslim", "mixtools")){
  missing <- setdiff(requirements,
                     rownames(installed.packages()))
  
  if (length(missing) != 0) {
    install.packages(missing, dependencies = TRUE)
  }
  invisible(lapply(requirements, require, character.only = TRUE))
}
#-------------------------------------------------------------------------------
