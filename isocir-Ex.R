pkgname <- "isocir"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
library('isocir')

assign(".oldSearch", search(), pos = 'CheckExEnv')
cleanEx()
nameEx("CIRE")
### * CIRE

flush(stderr()); flush(stdout())

### Name: CIRE
### Title: Circular Isotonic Regresssion Estimator
### Aliases: CIRE
### Keywords: circular isotonic order CIRE

### ** Examples


# We consider the following data from the package that are random circular data:
data(cirdata)
circular_ordered_estimator <- CIRE(cirdata)
# We can take the vector of the CIRE estimators:
circular_ordered_estimator $CIRE
# And the SCE:
circular_ordered_estimator $SCE

# Random data with a more complex order:
CIRE(cirdata, groups=c(1,1,2,3,5,3,4,6))




cleanEx()
nameEx("cirdata")
### * cirdata

flush(stderr()); flush(stdout())

### Name: cirdata
### Title: Random Circular Data.
### Aliases: cirdata
### Keywords: datasets circular

### ** Examples

data(cirdata)



cleanEx()
nameEx("cirgenes")
### * cirgenes

flush(stderr()); flush(stdout())

### Name: cirgenes
### Title: A set of angular measurements from cell-cycle experiments with
###   genes.
### Aliases: cirgenes
### Keywords: datasets circular genes

### ** Examples

data(cirgenes)



cleanEx()
nameEx("cond.test")
### * cond.test

flush(stderr()); flush(stdout())

### Name: cond.test
### Title: Conditional Test for Contrasting Circular Order
### Aliases: cond.test
### Keywords: circular isotonic CIRE test

### ** Examples

data(cirdata)
# Example without replications and a general isotropic order:
cond.test(cirdata, groups=c(1,2,1,3,3,4,5,6), kappa=0.2)
# Example with replications and the isotropic order (by default):
data(datareplic)
cond.test(data=datareplic)



cleanEx()
nameEx("datareplic")
### * datareplic

flush(stderr()); flush(stdout())

### Name: datareplic
### Title: Random Circular Data with Replications.
### Aliases: datareplic
### Keywords: datasets circular

### ** Examples

data(datareplic)



cleanEx()
nameEx("isocir")
### * isocir

flush(stderr()); flush(stdout())

### Name: isocir
### Title: S3 Objects of Class isocir.
### Aliases: isocir is.isocir print.isocir
### Keywords: isocir

### ** Examples

data(cirdata)
x <- CIRE(cirdata)
print(x)
is.isocir(x)
class(x)

plot(x)
class(x)

# If you want to use the CIRE in other calculations you can obtain it as a vector:
unlist(x$CIRE)
# But be careful because this unclass and lost attributes! 

# To create a new object of class CIRE:
y <- isocir()



cleanEx()
nameEx("mrl")
### * mrl

flush(stderr()); flush(stdout())

### Name: mrl
### Title: Mean Resultant Length
### Aliases: mrl
### Keywords: circular resultant length CIRE

### ** Examples

data(datareplic)
mrl(datareplic)



cleanEx()
nameEx("plot.isocir")
### * plot.isocir

flush(stderr()); flush(stdout())

### Name: plot.isocir
### Title: S3 Method to Plot S3 Objects of Class isocir
### Aliases: plot.isocir
### Keywords: circular plot CIRE

### ** Examples


data(cirdata)
result<-CIRE(cirdata)
plot(result)
plot(result,option="cirmeans")



cleanEx()
nameEx("sce")
### * sce

flush(stderr()); flush(stdout())

### Name: sce
### Title: Sum of Circular Error
### Aliases: sce
### Keywords: sce circular error

### ** Examples

data(cirdata)
exampledata1 <- cirdata
exampledata2 <- (cirdata+(pi/4))
sce(exampledata1,exampledata2)



### * <FOOTER>
###
cat("Time elapsed: ", proc.time() - get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
