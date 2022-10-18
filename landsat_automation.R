library(reticulate)
py_config()

Sys.setenv(RETICULATE_PYTHON = '/usr/local/lib/python2.7/')
#use_python('/usr/local/lib/python2.7/')

repl_python()

# ~~~~ Begin Python environment


import ee

#automate
folder <- 'LE07_L2SP_012054_20120217_20200909_02_T1'
bands <- c('_SR_B1') #etc
paste(folder, bands)

x <- list(r1, r2)
names(x) <- c("x", "y")
x$filename <- 'test.tif'
x$overwrite <- TRUE
m <- do.call(merge, x)











#how to deal if the UTM coordinates span zones?

labels <- read.csv('./data/training.csv' )
training <- getKMLcoordinates('./data/training.kml', ignoreAltitude=TRUE)
names(training) <- labels$Name
training <- training[5:length(training)]

training <- lapply(training, as.data.frame)
for (i in 1:length(training)){
  training[[i]]$label <- names(training[i])
}
training <- do.call(rbind, training)
colnames(training) <- c('x','y', 'label')

utm <- LongLatToUTM(training$x,training$y,17)
training <- cbind(utm, training)

training <- training %>%
  rename(x = 1,
         y = 2,
         long = 3,
         lat = 4)
training$label = as.factor(tolower(training$label))

training <- training[c(5:nrow(training)), -c(5,7,8,9)]
saveRDS(training, 'training.rds')










