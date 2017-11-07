### Datras download code

getDATRAS <- function(record, survey, startyear, endyear, quarters,
                      parallel = FALSE, cores = NULL, keepTime = FALSE) {

  if(keepTime == TRUE) strt <- Sys.time()
  #
  seqYear <- startyear:endyear
  #
  if(!record %in% c("HL", "HH", "CA")) stop("Please specify record type:
                                            HH (haul meta-data),
                                            HL (Species length-based data),
                                            CA (species age-based data)")
  getURL <- apply(expand.grid(record, survey, seqYear, quarters),
                  1,
                  function(x) paste0("http://datras.ices.dk/WebServices/",
                                     "DATRASWebService.asmx/get", x[1],
                                     "data?survey=", x[2],
                                     "&year=", x[3],
                                     "&quarter=", x[4]))
  #
  if(parallel == TRUE) {
    if(missing("cores")) stop("Please specify the number of cores.")
    #
    unregister <- function() {
      env <- foreach:::.foreachGlobals
      rm(list=ls(name=env), pos=env)
    } # close unregister
    #
    cl <- makeCluster(cores)
    registerDoParallel(cores = cl)
    
    #
    getDATA <- foreach(temp = getURL,
                       .combine = function(...) rbindlist(list(...), fill = TRUE),
                       .multicombine = T,
                       .inorder = F,
                       .maxcombine = 1000,
                       .packages = c("XML", "data.table")) %dopar% {
                         data.table(t(xmlSApply(xmlRoot(xmlTreeParse(temp, isURL = T,
                                                                     options = HUGE,
                                                                     useInternalNodes = T)),
                                                function(x) xmlSApply(x, xmlValue))))
                       } # close foreach %dopar%
    stopCluster(cl)
    unregister()
  } # close parallel == TRUE
  #
  if(parallel == FALSE) {
    getDATA <- foreach(temp = getURL,
                       .combine = function(...) rbindlist(list(...), fill = TRUE),
                       .multicombine=T,
                       .inorder=F,
                       .maxcombine=1000,
                       .packages = c("XML", "data.table")) %do% {
                         data.table(t(xmlSApply(xmlRoot(xmlTreeParse(temp, isURL = T,
                                                                     options = HUGE,
                                                                     useInternalNodes = T)),
                                                function(x) xmlSApply(x, xmlValue))))
                       } # close foreach %do%
  } # close parallel == FALSE
  if(keepTime == TRUE) print(Sys.time()-strt)
  return(getDATA)
} # close function

# Next, load the following packages.
library(XML, quietly = TRUE)
library(doParallel, quietly = TRUE)
library(parallel, quietly = TRUE)
library(foreach, quietly = TRUE)
library(data.table, quietly = TRUE)

haulData <- getDATRAS(record = "HH",
                      survey="NS-IBTS",
                      startyear = 2010,
                      endyear = 2011,
                      quarters = c(1,3),
                      parallel = TRUE,
                      cores = 4)
dim(haulData)


lengthData <- getDATRAS(record = "HL",
                        survey="BITS",
                        startyear = 2000,
                        endyear = 2001,
                        quarters = 1,
                        parallel = TRUE,
                        cores = 4)
dim(lengthData)

#Age Data
#To download the age-based information from all quarters in 1999 of the French Southern Atlantic
#Bottom Trawl Survey, use the following code:
  ageData <- getDATRAS(record = "CA",
                       survey="EVHOE",
                       startyear = 1999,
                       endyear = 1999,
                       quarters = c(1:4),
                       parallel = TRUE,
                       cores = 4)
dim(ageData)

