### * <HEADER>
###
attach(NULL, name = "CheckExEnv")
assign(".CheckExEnv", as.environment(2), pos = length(search())) # base
## add some hooks to label plot pages for base and grid graphics
setHook("plot.new", ".newplot.hook")
setHook("persp", ".newplot.hook")
setHook("grid.newpage", ".gridplot.hook")

assign("cleanEx",
       function(env = .GlobalEnv) {
	   rm(list = ls(envir = env, all.names = TRUE), envir = env)
           RNGkind("default", "default")
	   set.seed(1)
   	   options(warn = 1)
	   assign("T", delay(stop("T used instead of TRUE")),
		  pos = .CheckExEnv)
	   assign("F", delay(stop("F used instead of FALSE")),
		  pos = .CheckExEnv)
	   sch <- search()
	   newitems <- sch[! sch %in% .oldSearch]
	   for(item in rev(newitems))
               eval(substitute(detach(item), list(item=item)))
	   missitems <- .oldSearch[! .oldSearch %in% sch]
	   if(length(missitems))
	       warning("items ", paste(missitems, collapse=", "),
		       " have been removed from the search path")
       },
       env = .CheckExEnv)
assign("..nameEx", "__{must remake R-ex/*.R}__", env = .CheckExEnv) # for now
assign("ptime", proc.time(), env = .CheckExEnv)
grDevices::postscript("cgh-Examples.ps")
assign("par.postscript", graphics::par(no.readonly = TRUE), env = .CheckExEnv)
options(contrasts = c(unordered = "contr.treatment", ordered = "contr.poly"), pager="console")
library('cgh')

assign(".oldSearch", search(), env = .CheckExEnv)
assign(".oldNS", loadedNamespaces(), env = .CheckExEnv)
cleanEx(); ..nameEx <- "sw"

### * sw

flush(stderr()); flush(stdout())

### Name: sw
### Title: Perform the Smith-Waterman Algorithm
### Aliases: sw
### Keywords: misc

### ** Examples

## simluate vector of logratios
set.seed(3)
logratio <- c(rnorm(20) - 1, rnorm(20))

## invert sign of values and subtract threshold to ensure negative mean
x <- sw.threshold(logratio, function(x) median(x) + .2 * mad(x), sign = -1)

## perform Smith-Waterman algorithm
sw(x, trace = TRUE)
  


cleanEx(); ..nameEx <- "sw.perm.test"

### * sw.perm.test

flush(stderr()); flush(stdout())

### Name: sw.perm.test
### Title: Permutation Test for Smith-Waterman Algorithm
### Aliases: sw.perm.test
### Keywords: misc

### ** Examples

## simluate vector of logratios
set.seed(3)
logratio <- c(rnorm(20) - 1, rnorm(20))

## invert sign of values and subtract threshold to ensure negative mean
x <- sw.threshold(logratio, function(x) median(x) + .2 * mad(x), -1)

## perform Smith-Waterman
sw(x)

## perform permutation test on the islands identified
sw.perm.test(x, max.nIslands = NULL, nIter= 1e4)
  


cleanEx(); ..nameEx <- "sw.plot"

### * sw.plot

flush(stderr()); flush(stdout())

### Name: sw.plot
### Title: Plot Results of Smith-Waterman Algorithm
### Aliases: sw.plot
### Keywords: misc

### ** Examples

## simluate vector of logratios
set.seed(3)
logratio <- c(rnorm(20) - 1, rnorm(20))

## invert sign of values and subtract threshold to ensure negative mean
x <- sw.threshold(logratio, function(x) median(x) + .2 * mad(x), -1)

## perform permuation test for islands identified
p <- sw.perm.test(x, max.nIslands = NULL, nIter = 1e4)

## calculate robustness scores
r <- sw.rob(x)

## plot results
sw.plot(logratio, seq(length(logratio)),
  function(x) median(x) + .2 * mad(x), sign = -1, rob = r,
  main = paste("Toy dataset, highest-scoring island p =", p[1]))
  


cleanEx(); ..nameEx <- "sw.rob"

### * sw.rob

flush(stderr()); flush(stdout())

### Name: sw.rob
### Title: Robustness Calculation for Smith-Waterman Algorithm
### Aliases: sw.rob
### Keywords: misc

### ** Examples


## simluate vector of logratios
set.seed(3)
logratio <- c(rnorm(20) - 1, rnorm(20))

## invert sign of values and subtract threshold to ensure negative mean
x <- sw.threshold(logratio, function(x) median(x) + .2 * mad(x), -1)

## calculate robustness values
sw.rob(x)
  


cleanEx(); ..nameEx <- "sw.threshold"

### * sw.threshold

flush(stderr()); flush(stdout())

### Name: sw.threshold
### Title: Threshold function
### Aliases: sw.threshold
### Keywords: misc

### ** Examples

## simluate vector of logratios
set.seed(3)
logratio <- c(rnorm(20) - 1, rnorm(20))

## invert sign of values and subtract threshold to ensure negative mean
x <- sw.threshold(logratio, function(x) median(x) + .2 * mad(x), sign = -1)

## perform Smith-Waterman algorithm
sw(x, trace = TRUE)
  


### * <FOOTER>
###
cat("Time elapsed: ", proc.time() - get("ptime", env = .CheckExEnv),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
