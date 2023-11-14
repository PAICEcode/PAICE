pkgname <- "PAICE"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
library('PAICE')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("CmonsData")
### * CmonsData

flush(stderr()); flush(stdout())

### Name: CmonsData
### Title: Occurrence matrix of _Cistus monspeliensis_ in the Canary
###   Islands
### Aliases: CmonsData
### Keywords: datasets

### ** Examples

data(CmonsData)
CmonsData # Show data frame



cleanEx()
nameEx("CmonsNetwork")
### * CmonsNetwork

flush(stderr()); flush(stdout())

### Name: CmonsNetwork
### Title: Genealogical relationship of _Cistus monspeliensis_ haplotypes
### Aliases: CmonsNetwork
### Keywords: datasets

### ** Examples

data(CmonsNetwork)
CmonsNetwork # Show data frame



cleanEx()
nameEx("CmonsRare")
### * CmonsRare

flush(stderr()); flush(stdout())

### Name: CmonsRare
### Title: Simulated rarefaction curves of _Cistus monspeliensis_
### Aliases: CmonsRare
### Keywords: datasets

### ** Examples

data(CmonsRare)
str(CmonsRare) # Structure of data



cleanEx()
nameEx("PAICE")
### * PAICE

flush(stderr()); flush(stdout())

### Name: PAICE-package
### Title: Phylogeographic Analysis of Island Colonization Events
### Aliases: PAICE PAICE-package
### Keywords: package

### ** Examples

# Inference of minimum number of inter-island colonization events
data(CmonsData)
data(CmonsNetwork)
col <- colonization(data = CmonsData, network = CmonsNetwork)
col
summary(col)

# Asumptotic estimators of colonization events
# 25 replicates used in each sampling variable
## Note: The code is commented because 'CmonsRare' exists as an example
#set.seed(31)
#CmonsRare <- rarecol(data = CmonsData, network = CmonsNetwork,
#    replicates_field = 25, replicates_genetic = 25, monitor = TRUE,
#    mode = c(1, 2))
maxcol <- maxCol(data = CmonsRare, level = 0.95, del = 0.5, method = 1)
maxcol
summary(maxcol)

# Plotting results
par(mfrow = c(2, 2))
plot(CmonsRare)
par(fig = c(0, 1, 0, 0.5), new = TRUE)
plot(maxcol)



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("colonizations")
### * colonizations

flush(stderr()); flush(stdout())

### Name: colonization
### Title: Inference of minimum number of colonization events
### Aliases: colonization
### Keywords: univar

### ** Examples

data(CmonsData)
data(CmonsNetwork)
col <- colonization(data = CmonsData, network = CmonsNetwork)
col # Total of colonization events inferred
summary(col) # Detailed description of inferred colonization events



cleanEx()
nameEx("geneticResampling")
### * geneticResampling

flush(stderr()); flush(stdout())

### Name: geneticResampling
### Title: Simulate genetic sampling effort reduction
### Aliases: geneticResampling
### Keywords: datagen

### ** Examples

data(CmonsData)
data(CmonsNetwork)
# Delete position 462 of Cistus monspeliensis data
newdata <- geneticResampling(CmonsData, CmonsNetwork, 462)
newdata$data # New presences matrix of haplotypes
newdata$network # New genealogy



cleanEx()
nameEx("maxCol")
### * maxCol

flush(stderr()); flush(stdout())

### Name: maxCol
### Title: Asymptotic estimation of the number of colonization events
### Aliases: maxCol
### Keywords: univar

### ** Examples

# Use 'CmonsRare' data, a dataset generated using 25 replicates
# in both genetic and field sampling
data(CmonsRare)
maxcol <- maxCol(data = CmonsRare)
maxcol # Number of colonization estimated in each curve
summary(maxcol) # Description of curves
plot(maxcol) # Plotting estimations
# Plot all the information
par(mfrow = c(2, 2))
plot(CmonsRare) # First two plots with rarefaction curves
par(fig = c(0, 1, 0, 0.5), new = TRUE)
plot(maxcol) # Third plot with estimations



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("plot.maxCol")
### * plot.maxCol

flush(stderr()); flush(stdout())

### Name: plot.maxCol
### Title: Plot asymptotic estimators of colonization events
### Aliases: plot.maxCol
### Keywords: aplot

### ** Examples

# Use 'CmonsRare' data, a dataset generated using 25 replicates
# in both genetic and field sampling
data(CmonsRare)
maxcol <- maxCol(data = CmonsRare)
plot(maxcol)



cleanEx()
nameEx("plot.rarecol")
### * plot.rarecol

flush(stderr()); flush(stdout())

### Name: plot.rarecol
### Title: Plot rarefaction curve of colonization events
### Aliases: plot.rarecol
### Keywords: aplot

### ** Examples

# Use 'CmonsRare' data, a dataset generated using 25 replicates
# in both genetic and field sampling
data(CmonsRare)
plot(CmonsRare)



cleanEx()
nameEx("rarecol")
### * rarecol

flush(stderr()); flush(stdout())

### Name: rarecol
### Title: Rarefaction curve of colonization events
### Aliases: rarecol
### Keywords: datagen

### ** Examples

data(CmonsData)
data(CmonsNetwork)
# Build rarefaction curves with 5 field and genetic replicates
## Note: more replicates are needed to build accurate curves
## Note: 5 replicates are relatively fast and adequate to
##       explore the data
rcol <- rarecol(data = CmonsData, network = CmonsNetwork,
                replicates_field = 5, replicates_genetic = 5,
                monitor = TRUE, mode = c(1, 2))
par(mfrow = c(1, 2))
plot(rcol) # Plotting results



graphics::par(get("par.postscript", pos = 'CheckExEnv'))
cleanEx()
nameEx("read.rarecol")
### * read.rarecol

flush(stderr()); flush(stdout())

### Name: read.rarecol
### Title: Read files containing rarefaction curves of colonization events
### Aliases: read.rarecol
### Keywords: file

### ** Examples

data(CmonsData)
data(CmonsNetwork)
# Make rarefaction curves and save it in working directory,
## Note: only one replicate per sampling to it quickly
rarecol(data = CmonsData, network = CmonsNetwork,
        replicates_field = 1, replicates_genetic = 1,
        monitor = TRUE, file = "rareData")
# Genetic estimation has the suffix "_gen" and the field "_field"
raredata <- read.rarecol(gen = "rareData_gen.csv",
                         field = "rareData_field.csv")
str(raredata) # Show structure of data imported
# Remove files created
file.remove("rareData_gen.csv", "rareData_field.csv")



### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
