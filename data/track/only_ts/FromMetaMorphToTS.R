setwd("~/Dropbox/CHUKNORRIS/data/track/only_ts/")
fls <- list.files("~/Dropbox/CHUKNORRIS/data/track/")
for(i in c(1:3,5:7)){
  dat <- read.csv(paste("~/Dropbox/CHUKNORRIS/data/track/", fls[i], sep = ""))
  for(j in 1:3){
    ts.save <- cbind(cumsum(dat$Time.Interval)/60,dat$Velocity[which(dat$Object..== j)]*60*0.2160014)
    colnames(ts.save) <- c("Time","Growth rate")
    write.csv(ts.save,paste("TS_Object",j,fls[i],sep="_"), row.names = FALSE)
    
  }
  
}

i <- 4
dat <- read.csv(paste("~/Dropbox/CHUKNORRIS/data/track/", fls[i], sep = ""))
ys <- (60 * 0.2160014 * dat$X..Velocity. * dat$X..Time.Interval.) / 4
xs <- cumsum(c(0, rep(4, dim(dat)[1] - 1))) / 60
ts.save <- cbind(xs, ys)
colnames(ts.save) <- c("Time","Growth rate")
write.csv(ts.save,paste("TS_Object",1,fls[i],sep="_"), row.names = FALSE)
