estDailyFromSurfaces <- function(eList, localsurfaces = NA, localDaily = NA) {
  ## localsurfaces = NA
  ## localDaily = NA

  if(!is.egret(eList)){
    stop("Please check eList argument")
  }

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
  allLogQsReplicated <- allLogQsByDayOfYear[as.character(localDaily$Day)]

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
  allLogQsByDayOfYear <- split(localDaily$LogQ, localDaily$Day)

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
                       Sample = eList$Sample,
                       windowY=windowY, windowQ=windowQ, windowS=windowS,
                       minNumObs=minNumObs, minNumUncen=minNumUncen,edgeAdjust=edgeAdjust,
                       verbose=verbose)

  eList$Sample <- Sample1

  if(verbose) cat("\nNext step running  estSurfaces with survival regression:\n")

  surfaces1 <- estSurfaces(eList,
                         windowY = windowY, windowQ=windowQ, windowS=windowS,
                         minNumObs = minNumObs, minNumUncen = minNumUncen, edgeAdjust = edgeAdjust,
                         verbose = verbose, run.parallel = run.parallel)

  eList$surfaces <- surfaces1

  Daily1 <- estDailyFromSurfaces(eList)

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

estCrossVal<-function(DecLow,DecHigh, Sample, windowY = 7, windowQ = 2,
                      windowS = 0.5, minNumObs = 100, minNumUncen = 50,
                      edgeAdjust=TRUE, verbose = TRUE){
  #  this function fits the WRTDS model making an estimate of concentration for every day
  #    But, it uses leave-one-out-cross-validation
  #    That is, for the day it is estimating, it leaves that observation out of the data set
  #      It returns a Sample data frame with three added columns
  #      yHat, SE, and ConcHat

  localSample <- Sample
  originalColumns <- names(localSample)
  numObs<-nrow(localSample)
  yHat<-rep(0,numObs)
  SE<-rep(0,numObs)
  ConcHat<-rep(0,numObs)
  iCounter<-seq(1,numObs)

  if(verbose) cat("\n estCrossVal % complete:\n")

  colToKeep <- c("ConcLow","ConcHigh","Uncen","DecYear","SinDY","CosDY","LogQ")
  SampleCrossV <- localSample[,which(originalColumns %in% colToKeep)]

  SampleCV<-data.frame(SampleCrossV,iCounter,yHat,SE,ConcHat)

  printUpdate <- floor(seq(1,numObs,numObs/100))
  endOfLine <- seq(10,100,10)

  for(i in 1:numObs) {
    if(i %in% printUpdate & verbose) {
      cat(floor(i*100/numObs),"\t")
      if (floor(i*100/numObs) %in% endOfLine) cat("\n")
    }

    SampleMinusOne<-SampleCV[SampleCV$iCounter!=i,]

    result <- runSurvReg(estPtYear = SampleCrossV$DecYear[i],
                       estPtLQ = SampleCrossV$LogQ[i],
                       DecLow = DecLow,
                       DecHigh = DecHigh,
                       Sample = SampleMinusOne,
                       windowY,windowQ,windowS,
                       minNumObs,minNumUncen,
                       edgeAdjust=edgeAdjust, verbose=FALSE, run.parallel = FALSE)

    yHat[i]<-result[1]
    SE[i]<-result[2]
    ConcHat[i]<-result[3]

  }

  localSample$yHat <- yHat
  localSample$SE <- SE
  localSample$ConcHat <- ConcHat
  SampleCrossV <- localSample
  return(SampleCrossV)
}

runSurvReg<-function(estPtYear,estPtLQ,DecLow,DecHigh,Sample,
                     windowY=7, windowQ=2, windowS=0.5,
                     minNumObs=100, minNumUncen=50, verbose = TRUE,interactive=NULL,
                     edgeAdjust=TRUE, run.parallel = FALSE) {

  if(!is.null(interactive)) {
    warning("The argument 'interactive' is deprecated. Please use 'verbose' instead")
    verbose <- interactive
  }

  localSample <- Sample
  if(any(is.na(localSample$LogQ))){
    message("Removing Sample data that does not have corresponding flow data")
    localSample <- localSample[!is.na(localSample$LogQ),]
  }
  numSamples <- length(localSample$DecYear)

  numEstPt<-length(estPtYear)

  printUpdate <- floor(seq(1,numEstPt,numEstPt/100))
  endOfLine <- seq(10,100,10)

  if (minNumUncen >= sum(localSample$Uncen)) stop('minNumUncen is greater than total number of samples')
  if (minNumObs >= nrow(localSample)) stop('minNumObs is greater than total number of samples')

  warningFlag <- 0
  n <- NULL

  if(run.parallel){
    `%dopar%` <- foreach::`%dopar%`
    wrtds_return_list <- foreach::foreach(n = 1:numEstPt, .packages=c('EGRET')) %dopar% {
                      wrtd27] -10.481699  -9.419595  -8.357491  -7.295387  -6.233283  -5.171179
[433]  -4.109075  -3.046971 -16.854322 -15.792218 -14.730115 -13.668011
[439] -12.605907 -11.543803 -10.481699  -9.419595  -8.357491  -7.295387
[445]  -6.233283  -5.171179  -4.109075  -3.046971 -16.854322 -15.792218
[451] -14.730115 -13.668011 -12.605907 -11.543803 -10.481699  -9.419595
[457]  -8.357491  -7.295387  -6.233283  -5.171179  -4.109075  -3.046971
> s_returns <- run_WRTDS(estY = estPtYear[n], estLQ = estPtLQ[n],
                                                 localSample =localSample,DecLow = DecLow,DecHigh = DecHigh,
                                                 minNumObs = minNumObs,minNumUncen = minNumUncen,
                                                 windowY = windowY, windowQ = windowQ, windowS = windowS,
                                                 edgeAdjust = edgeAdjust)
                    }

    warningFlag <- sum(sapply(wrtds_return_list, function(x) x[["warningFlag"]]))
    resultSurvReg <- t(sapply(wrtds_return_list, function(x) x[["survReg"]]))

  } else {
    resultSurvReg<-array(0,c(numEstPt,3))
    if (verbose) cat("Survival regression (% complete):\n")

    for (i in 1:numEstPt) {
      print(i)
      print('printing estPtLQ[i]')
      print(estPtLQ[i])

      wrtds_return <- run_WRTDS(estY = estPtYear[i],
                                estLQ = estPtLQ[i],
                                localSample = localSample,
                                DecLow = DecLow,DecHigh = DecHigh,
                                minNumObs = minNumObs,minNumUncen = minNumUncen,
                                windowY = windowY, windowQ = windowQ, windowS = windowS,
                                edgeAdjust = edgeAdjust)

      if (i %in% printUpdate & verbose) {
        cat(floor(i*100/numEstPt),"\t")
        if (floor(i*100/numEstPt) %in% endOfLine) cat("\n")
      }
      warningFlag <- warningFlag + wrtds_return$warningFlag
      resultSurvReg[i,] <- wrtds_return$survReg
    }
  }

  if (warningFlag > 0){
    message("\nIn model estimation, the survival regression function was run ", numEstPt, " times (for different combinations of discharge and time).  In ", warningFlag, " of these runs it did not properly converge. This does not mean that the model is unacceptable, but it is a suggestion that there may be something odd about the data set. You may want to check for outliers, repeated values on a single date, or something else unusual about the data.")
  }

  if (verbose) cat("\nSurvival regression: Done")

  return(resultSurvReg)
}

run_WRTDS <- function(estY, estLQ,
                      localSample,DecLow,DecHigh,
                      minNumObs,minNumUncen,
                      windowY, windowQ, windowS,
                      edgeAdjust){
  # This loop takes us through all the estimation points
  # We always reset the window widths to their starting values, because they may
  #   have been widened in the process
  tempWindowY<-windowY
  tempWindowQ<-windowQ
  tempWindowS<-windowS

  distLow <- estY-DecLow
  distHigh <- DecHigh-estY

  survReg <- c(NA, NA, NA)
  warningFlag <- 0

  if(all(is.na(c(distLow,distHigh)))){
    return(list(survReg=survReg, warningFlag=warningFlag))
  }

  distTime <- min(distLow,distHigh)

  if (edgeAdjust & !is.na(distTime)) {
    tempWindowY <- if(distTime>tempWindowY) tempWindowY else ((2 * tempWindowY) - distTime)
  }

  k <- 1

  repeat{
    #  We subset the sample frame by time, to narrow the set of data to run through in the following steps

    Sam <- localSample[abs(localSample$DecYear-estY) <= tempWindowY,]
    diffY<-abs(Sam$DecYear-estY)
    weightY<-triCube(diffY,tempWindowY)
    weightQ<-triCube(Sam$LogQ-estLQ,tempWindowQ)
    diffUpper<-ceiling(diffY)
    diffLower<-floor(diffY)
    diffSeason<-pmin(abs(diffUpper-diffY),abs(diffY-diffLower))
    weightS<-triCube(diffSeason,tempWindowS)
    Sam$weight<-weightY*weightQ*weightS
    Sam<-subset(Sam,weight>0)
    numPosWt<-length(Sam$weight)
    numUncen<-sum(Sam$Uncen)
    tempWindowY<-tempWindowY*1.1
    tempWindowQ<-tempWindowQ*1.1
    k <- k + 1
    if(k > 10000) message("Problems converging")
    # the next line is designed so that if the user sets windowS so it includes
    # data from all seasons, the widening process leaves it alone
    tempWindowS<-if(windowS<=0.5) min(tempWindowS*1.1,0.5) else windowS
    if(numPosWt>=minNumObs&numUncen>=minNumUncen | k > 10000) break
  }

  # now we are ready to run Survival Regression
  weight<-Sam$weight
  aveWeight<-sum(weight)/numPosWt
  weight<-weight/aveWeight
  Sam <- data.frame(Sam)

  x <- tryCatch({
    survModel <- survival::survreg(survival::Surv(log(ConcLow),log(ConcHigh),type="interval2") ~
                                     DecYear+LogQ+SinDY+CosDY,data=Sam,weights=weight,dist="gaus")
  }, warning=function(w) {

    if(w$message == "Ran out of iterations and did not converge"){
      Sam2 <- jitterSam(Sam)
    } else {
      Sam2 <- Sam
    }

    survModel <- survival::survreg(survival::Surv(log(ConcLow),log(ConcHigh),type="interval2") ~
                                     DecYear+LogQ+SinDY+CosDY,data=Sam2,weights=weight,dist="gaus")

    return(survModel)
  }, error=function(e) {
    message(e, "Error")
    return(NULL)
  })

  if(class(x) == "survreg") {
    newdf<-data.frame(DecYear=estY,LogQ=estLQ,SinDY=sin(2*pi*estY),CosDY=cos(2*pi*estY))
    #   extract results at estimation point
    yHat<-predict(x,newdf)
    SE<-x$scale
    bias<-exp((SE^2)/2)
    survReg[1]<-yHat
    survReg[2]<-SE
    survReg[3]<-bias*exp(yHat)
  }

  if(all(is.na(x))){
    warningFlag <- 1
  }
  return(list(survReg=survReg, warningFlag=warningFlag))
}

estSurfaces<-function(eList, surfaceStart=NA, surfaceEnd=NA, localSample=NA,
                      windowY=7,windowQ=2,windowS=0.5,
                      minNumObs=100,minNumUncen=50,edgeAdjust=TRUE,
                      verbose = TRUE, interactive=NULL,
                      run.parallel = FALSE){

  # this function estimates the 3 surfaces based on the Sample data
  # one is the estimated log concentration (yHat)
  # the second is the estimated standard error (SE)
  # the third is the estimated concentration (ConcHat)
  # they are mapped as an array that covers the complete space of daily discharge and time
  # the first index is discharge, layed out in 14 equally spaced levels of log(Q)
  # the second index is time, layed out as 16 increments of the calendar year, starting January 1.
  # it returns the data frame called surfaces
  #
  if(!is.null(interactive)) {
    warning("The argument 'interactive' is deprecated. Please use 'verbose' instead")
    verbose <- interactive
  }
  if(!is.egret(eList)){
    stop("Please check eList argument")
  }
  localINFO <- getInfo(eList)
  localDaily <- getDaily(eList)

  if(all(is.na(localSample))){
    localSample <- eList$Sample
  }

  highLow <- decimalHighLow(localSample)

  DecHigh <- highLow[["DecHigh"]]
  DecLow <- highLow[["DecLow"]]

  surfaceInfo <- surfaceIndex(localDaily)
  vectorYear <- surfaceInfo[['vectorYear']]
  vectorLogQ <- surfaceInfo[['vectorLogQ']]

  LogQ <- seq(surfaceInfo[['bottomLogQ']], by=surfaceInfo[['stepLogQ']], length.out=surfaceInfo[['nVectorLogQ']])

  if(is.na(surfaceStart) && is.na(surfaceEnd)){

    nVectorYear<-length(vectorYear)
    estPtYear<-rep(vectorYear, each=14)

    Year <- seq(surfaceInfo[['bottomYear']], by=surfaceInfo[['stepYear']], length.out=surfaceInfo[['nVectorYear']])

  } else {

    sliceIndex <- which(vectorYear >= decimalDate(as.Date(surfaceStart)) & vectorYear <=
                          decimalDate(as.Date(surfaceEnd)))
    Year <- vectorYear[c(sliceIndex[1]-1, sliceIndex, tail(sliceIndex, n = 1)+1)]

    nVectorYear <- length(Year)
    estPtYear <- rep(Year,each=14)

  }

  estPtLogQ<-rep(vectorLogQ,nVectorYear)

  resultSurvReg<-runSurvReg(estPtYear, estPtLogQ,
                            DecLow, DecHigh, localSample,
                            windowY,windowQ,windowS,
                            minNumObs,minNumUncen,
                            edgeAdjust=edgeAdjust,
                            verbose = verbose,run.parallel=run.parallel)

  surfaces<-array(0,dim=c(14,nVectorYear,3))

  for(iQ in 1:14) {
    for(iY in 1:nVectorYear){
      k<-(iY-1)*14+iQ
      surfaces[iQ,iY,]<-resultSurvReg[k,]
    }
  }

  attr(surfaces, "surfaceIndex") <- surfaceInfo
  attr(surfaces, "LogQ") <- LogQ
  attr(surfaces, "Year") <- Year

  return(surfaces)
}

surfaceIndex<-function(Daily){
  # this function contains the same code that comes at the start of
  # estSurfaces, it just computes the parameters of the grid
  # used for the surfaces so that they can be stored for future use
  # the first index is discharge, layed out in 14 equally spaced levels of log(Q)
  # the second index is time, layed out as 16 increments of the calendar year, starting January 1.
  #  Note: I don't think this is the smartest way to do this, but I'm not sure what to do here
  #  I don't like trying to have the same code twice
  #

  localDaily <- Daily

  bottomLogQ<- min(localDaily$LogQ, na.rm = TRUE) - 0.05
  topLogQ <- max(localDaily$LogQ, na.rm = TRUE) + 0.05
  stepLogQ <-(topLogQ-bottomLogQ)/13
  vectorLogQ <- seq(bottomLogQ,topLogQ,stepLogQ)
  stepYear<-1/16
  bottomYear<-floor(min(localDaily$DecYear, na.rm = TRUE))
  topYear<-ceiling(max(localDaily$DecYear, na.rm = TRUE))
  vectorYear<-seq(bottomYear,topYear,stepYear)
  nVectorYear<-length(vectorYear)

  surfaceIndexParameters<-list(bottomLogQ=bottomLogQ,
                            stepLogQ=stepLogQ,
                            nVectorLogQ=14,
                            bottomYear=bottomYear,
                            stepYear=stepYear,
                            nVectorYear=nVectorYear,
                            vectorYear=vectorYear,
                            vectorLogQ=vectorLogQ)
  return(surfaceIndexParameters)
}
