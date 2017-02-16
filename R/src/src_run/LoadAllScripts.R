#-------------------------------------------------------------------------------
LoadAllScripts <- function(src.path = "~/Dropbox/CHUKNORRIS/R/src/"){
  all.scripts <- list.files(src.path, pattern = "[.]R$", recursive = TRUE)
  invisible(lapply(paste(src.path, all.scripts, sep = ""), source))
}
#-------------------------------------------------------------------------------
