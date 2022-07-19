decimalDateWY <- function(dates, wy_type = 'usgs') {
    yd <- index(dates) - 1
    ## yend <- yday(as.Date(paste0(year(dates), "-12-31")))
    yend <- rep(length(dates), length(dates))

    ydec <- as.integer(as.character(water_year(dates, wy_type))) + yd/yend

  return(ydec)
}

estDailyFromSurfaces <- function(eList, localsurfaces = NA, localDaily = NA) {

  if(!is.egret(eList)){
    stop("Please check eList argument")
  }

  print('estDailyFromSurfaces')

  localDaily <- getSurfaceEstimates(eList, localsurfaces=localsurfaces, localDaily = localDaily)
  # Calculate "flow-normalized" concentration and flux:
  allLogQsByDayOfYear <- bin_Qs(localDaily)

  concFlux_list <- getConcFluxFromSurface(eList, allLogQsByDayOfYear, localDaily = localDaily, localsurfaces=localsurfaces)

  ## cFx <- concFlux_list[["allConcReplicated"]]
  ## dFx <- concFlux_list[["allDatesReplicated"]]
  ## test <- tapply(concFlux_list[["allConcReplicated"]], concFlux_list[["allDatesReplicated"]], "mean")

  # Finally bin the collective results by days (the decimal year), and calculate the desired means.
  localDaily$FNConc <-  as.numeric(tapply(concFlux_list[["allConcReplicated"]], concFlux_list[["allDatesReplicated"]], "mean"))
  localDaily$FNFlux <-  as.numeric(tapply(concFlux_list[["allFluxReplicated"]], concFlux_list[["allDatesReplicated"]], "mean"))

  return(localDaily)
}

getConcFluxFromSurface <- function(eList, allLogQsByDayOfYear, localDaily, localsurfaces = NA){

  if(all(is.na(localsurfaces))){
    localsurfaces <- getSurfaces(eList)
  }


  print('getConcFluxFromSurface')

  # First argument in calls below is the "known" x-y-z surface, second argument is matrix of
  # "target" x-y points.
  if("LogQ" %in% names(attributes(localsurfaces))){
    LogQ <- attr(localsurfaces, "LogQ")
  } else {
    localINFO <- getInfo(eList)
    if(all(c("bottomLogQ","stepLogQ","nVectorLogQ") %in% names(localINFO))){
      LogQ <- seq(localINFO$bottomLogQ, by=localINFO$stepLogQ, length.out=localINFO$nVectorLogQ)
    } else {
      surfaceIndexParameters<-surfaceIndex(eList$Daily)
      bottomLogQ<-surfaceIndexParameters[['bottomLogQ']]
      stepLogQ<-surfaceIndexParameters[['stepLogQ']]
      nVectorLogQ<-surfaceIndexParameters[['nVectorLogQ']]
      LogQ <- seq(bottomLogQ, by=stepLogQ, length.out=nVectorLogQ)
    }}

  if("Year" %in% names(attributes(localsurfaces))){
    Year <- attr(localsurfaces, "Year")
  } else {
    localINFO <- getInfo(eList)
    if(all(c("bottomYear","stepYear","nVectorYear") %in% names(localINFO))){
      Year <- seq(localINFO$bottomYear, by=localINFO$stepYear, length.out=localINFO$nVectorYear)
    } else {
      surfaceIndexParameters <- surfaceIndex(eList$Daily)
      bottomYear <- surfaceIndexParameters[['bottomYear']]
      stepYear <- surfaceIndexParameters[['stepYear']]
      nVectorYear <- surfaceIndexParameters[['nVectorYear']]
      Year <- seq(bottomYear, by=stepYear, length.out=nVectorYear)
    }}

  # Using the above data structure as a "look-up" table, list all LogQ values that occured on every
  # day of the entire daily record. When "unlisted" into a vector, these will become the "x" values
  # for the interpolation.
  allLogQsReplicated <- allLogQsByDayOfYear[index(localDaily$Day)]

  # Replicate the decimal year field for each day of the record to correspond to all the LogQ
  # values listed for that day. These are the "y" values for the interpolation.
  allDatesReplicated <- rep(localDaily$DecYear, lapply(allLogQsReplicated, length))

  # Interpolate.
  allConcReplicated <- fields::interp.surface( obj=list(x=LogQ,y=Year,z=localsurfaces[,,3]),
                                       loc=data.frame(	unlist(x=allLogQsReplicated),
                                                       y=allDatesReplicated))
  allFluxReplicated <- allConcReplicated * exp(as.numeric(unlist(allLogQsReplicated))) * 86.4

  return(list(allFluxReplicated=allFluxReplicated,
              allConcReplicated=allConcReplicated,
              allDatesReplicated=allDatesReplicated))

}

getSurfaceEstimates <- function(eList, localsurfaces=NA, localDaily = NA){

  if(all(is.na(localDaily))){
    localDaily <- getDaily(eList)
  }

  if(all(is.na(localsurfaces))){
    localsurfaces <- getSurfaces(eList)
  }

  print('getSurfaceEstimates')

  if("LogQ" %in% names(attributes(localsurfaces))){
    LogQ <- attr(localsurfaces, "LogQ")
  } else {
    localINFO <- getInfo(eList)
    LogQ <- seq(localINFO$bottomLogQ, by=localINFO$stepLogQ, length.out=localINFO$nVectorLogQ)
  }

  if("Year" %in% names(attributes(localsurfaces))){
    Year <- attr(localsurfaces, "Year")
  } else {
    localINFO <- getInfo(eList)
    Year <- seq(localINFO$bottomYear, by=localINFO$stepYear, length.out=localINFO$nVectorYear)
  }

  localDaily$yHat <- fields::interp.surface(obj=list(x=LogQ,y=Year,z=localsurfaces[,,1]),
                                    loc=data.frame(localDaily$LogQ, localDaily$DecYear))
  localDaily$SE <- fields::interp.surface(obj=list(x=LogQ,y=Year,z=localsurfaces[,,2]),
                                  loc=data.frame(localDaily$LogQ, localDaily$DecYear))
  localDaily$ConcDay <- fields::interp.surface(obj=list(x=LogQ,y=Year,z=localsurfaces[,,3]),
                                       loc=data.frame(localDaily$LogQ, localDaily$DecYear))
  localDaily$FluxDay <- as.numeric(localDaily$ConcDay * localDaily$Q * 86.4)

  return(localDaily)
}

bin_Qs <- function(localDaily){

  print('binQs')
  allLogQsByDayOfYear <- split(localDaily$LogQ, index(localDaily$Day))

  # account for leap day, day of year '60'
  ## allLogQsByDayOfYear[['59']] <- c(unlist(allLogQsByDayOfYear['59']),   # Bob's convention
  ##                                  unlist(allLogQsByDayOfYear['60']))
  ## allLogQsByDayOfYear['60'] <- allLogQsByDayOfYear['59']

  return(allLogQsByDayOfYear)

}

modelEstimation<-function(eList,
                          windowY = 7, windowQ = 2, windowS = 0.5,
                          minNumObs = 100, minNumUncen = 50,
                          edgeAdjust = TRUE, verbose = TRUE,
                          run.parallel = FALSE){

  print('----------------- MODEL ESTIMATION')

  if(!is.egret(eList)){
    stop("Please check eList argument")
  }

  eList <- setUpEstimation(eList = eList,
                           windowY = windowY, windowQ = windowQ, windowS = windowS,
                           minNumObs = minNumObs, minNumUncen = minNumUncen,
                           edgeAdjust = edgeAdjust, verbose = verbose)

  if(verbose) cat("\n first step running estCrossVal may take about 1 minute")

  Sample1 <- estCrossVal(eList$Daily$DecYear[1],
                       eList$Daily$DecYear[length(eList$Daily$DecYear)],
                       eList$Sample,
                       windowY=windowY, windowQ=windowQ, windowS=windowS,
                       minNumObs=minNumObs, minNumUncen=minNumUncen,edgeAdjust=edgeAdjust,
                       verbose=verbose)

  eList$Sample <- Sample1

  if(verbose) cat("\nNext step running  estSurfaces with survival regression:\n")

  surfaces1 <- estSurfaces(eList,
                         windowY = windowY, windowQ=windowQ, windowS=windowS,
                         minNumObs = minNumObs, minNumUncen = minNumUncen, edgeAdjust = edgeAdjust,
                         verbose = verbose, run.parallel = run.parallel)

  print('__ surface')
  eList$surfaces <- surfaces1

  Daily1 <- estDailyFromSurfaces(eList)

  print('__ daily')
  print(nrow(Daily1))
  eList$Daily <- Daily1

  checkSurfaceSpan(eList)

  return(eList)

}



setUpEstimation<-function(eList,
                          windowY=7, windowQ=2, windowS=0.5,
                          minNumObs=100,minNumUncen=50,
                          edgeAdjust=TRUE, verbose=TRUE, interactive = NULL){

  if(!is.null(interactive)) {
    warning("The argument 'interactive' is deprecated. Please use 'verbose' instead")
    verbose <- interactive
  }

  localINFO <- getInfo(eList)
  localSample <- getSample(eList)
  localDaily <- getDaily(eList)

  if(!all(c("Q","LogQ") %in% names(localSample))){
    eList <- mergeReport(INFO=localINFO, Daily = localDaily, Sample = localSample, verbose=verbose)
  }

  if(any(localSample$ConcLow[!is.na(localSample$ConcLow)] == 0)){
    stop("modelEstimation cannot be run with 0 values in ConcLow. An estimate of the reporting limit needs to be included. See fixSampleFrame to adjust the Sample data frame")
  }

  numDays <- length(localDaily$DecYear)
  DecLow <- localDaily$DecYear[1]
  DecHigh <- localDaily$DecYear[numDays]

  surfaceIndexParameters<-surfaceIndex(localDaily)
  localINFO$bottomLogQ<-surfaceIndexParameters[['bottomLogQ']]
  localINFO$stepLogQ<-surfaceIndexParameters[['stepLogQ']]
  localINFO$nVectorLogQ<-surfaceIndexParameters[['nVectorLogQ']]
  localINFO$bottomYear<-surfaceIndexParameters[['bottomYear']]
  localINFO$stepYear<-surfaceIndexParameters[['stepYear']]
  localINFO$nVectorYear<-surfaceIndexParameters[['nVectorYear']]
  localINFO$windowY<-windowY
  localINFO$windowQ<-windowQ
  localINFO$windowS<-windowS
  localINFO$minNumObs<-minNumObs
  localINFO$minNumUncen<-minNumUncen
  localINFO$numDays <- numDays
  localINFO$DecLow <- DecLow
  localINFO$DecHigh <- DecHigh
  localINFO$edgeAdjust <- edgeAdjust

  eList$INFO <- localINFO

  return(eList)

}
