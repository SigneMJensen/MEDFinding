amavg<-function (modelList, modelWeights, expressions, interval=c("Buckland", "Asymptotic"),level=0.95,
                 distribution="normal",seType="mod")
{
  require(sandwich, quietly = TRUE)
  require(car, quietly = TRUE)
  deltab<-function (object, g, func = g, ...)
  {
    if (!is.character(g))
      stop("The argument 'g' must be a character string")
    para <- coef(object)
    if(inherits(object,"lm")){
      para.names <- names(coef(object))
      para.names[1] <- gsub("\\(Intercept\\)", "Intercept", para.names[1])
      names(para) <- para.names
    }else{
      if(inherits(object,"drc")){
        coefVec<-coef(object)
        para.names<-sapply(strsplit(names(coefVec), ":"), "[[", 1)
        names(para)<-para.names
      }else{
        para.names <- names(para)
      }}
    g <- parse(text = g)
    q <- length(para)
    for (i in 1:q) {
      assign(names(para)[i], para[i])
    }
    gd <- rep(0,q)
    for (i in 1:q) {
      gd[i] <- eval(D(g, names(para)[i]))
    }
    gd
  }
  makeIIDdecomp <- function(modelObject,g)
  {
    numObsUsed <- ifelse(inherits(modelObject, "coxph"),
                         modelObject$n, ifelse(inherits(modelObject,"nls"),
                                               length(predict(modelObject)),ifelse(inherits(modelObject,"drc"),
                                                                                   length(predict(modelObject)), nrow(modelObject$model))))
    res.df<-df.residual(modelObject)
    db<-deltab(modelObject,g)
    iidVec0 <- db %*% bread(modelObject) %*% t(estfun(modelObject))
    moNAac <- modelObject$na.action
    numObs <- numObsUsed + length(moNAac)
    iidVec <- rep(0, numObs)
    if (!is.null(moNAac)) {
      iidVec[-moNAac] <- sqrt(numObs/numObsUsed) * iidVec0
    } else {
      iidVec <- iidVec0
    }
    list(iidVec = iidVec, numObsUsed = numObsUsed, numObs = numObs, res.df=res.df)
  }
  iidList <- list()
  numModels <- length(modelList)
  for(i in 1:numModels)
  {
    iidList[[i]]<- makeIIDdecomp(modelList[[i]], expressions[[i]])
  }
  iidresp <- matrix(as.vector(unlist(lapply(iidList, function(listElt) {
    listElt[[1]]}))), nrow = numModels, byrow = TRUE)
  numObsUsed <- as.vector(unlist(lapply(iidList, function(listElt) {
    listElt[[2]]})))
  res.df.vec <- as.vector(unlist(lapply(iidList, function(listElt) {
    listElt[[4]]})))
  thetaEst <- rep(NA, numModels)
  thetaSe <- rep(NA, numModels)
  for(i in 1:numModels)
  {
    if (inherits(modelList[[i]], "drc")){
      coefVec <- coef(modelList[[i]])
      names(coefVec) <- sapply(strsplit(names(coefVec), ":"), "[[", 1)
      deltaRes <- deltaMethod(coefVec,expressions[[i]],vcov(modelList[[i]]))
    } else {
      deltaRes <- deltaMethod(modelList[[i]],expressions[[i]])
    }
    thetaEst[i] <- deltaRes[1]
    thetaSe[i] <- deltaRes[2]
  }
  thetaEst <- unlist(thetaEst)
  thetaSe <- unlist(thetaSe)
  thetaMA <- as.numeric(modelWeights%*%thetaEst)
  if(identical(distribution, "t")){
    quant <- qt(1 - (1 - level)/2, df.residual(modelList[[1]]))
  } else {
    quant <- qnorm(1 - (1 - level)/2)
  }
  numObs <- iidList[[1]]$numObs
  covar <- (iidresp %*% t(iidresp)) / numObs
  vcMat <- covar / numObs # Defining the finite-sample variance-covariance matrix
  ## Replacing sandwich estimates by model-based standard errors
  modbas <- seType == "mod"
  if (any(modbas))
  {
    corMat <- cov2cor(vcMat)
    ## Retrieving standard errors for the specified estimate from the individual fits
    modSE <- thetaSe
    sanSE <- sqrt(diag(vcMat))
    sanSE[modbas] <- modSE[modbas]
    vcMat <- diag(sanSE,nrow=length(sanSE)) %*% corMat %*% diag(sanSE,nrow=length(sanSE))
  }
  if (identical(interval, "Asymptotic"))
  {
    numObs <- iidList[[1]]$numObs
    varMA <- modelWeights %*% t(modelWeights %*% vcMat)
    seMA <- sqrt(diag(varMA))
    quantVal <- quant * seMA
    retMat <- as.matrix(cbind(thetaMA, seMA, thetaMA - quantVal, thetaMA + quantVal))
    colnames(retMat) <- c("Estimate", "Std. Error", "Lower", "Upper")
  }
  if (identical(interval, "Buckland"))
  {
    seVec <- apply(sqrt(thetaSe^2 + (t(t(thetaEst) - mean(thetaMA)))^2) * modelWeights, 2,
                   sum)
    quantVal <- quant * seVec
    retMat <- as.matrix(cbind(thetaMA, seVec, thetaMA - quantVal, thetaMA + quantVal))
    colnames(retMat) <- c("Estimate", "Std. Error", "Lower", "Upper")
  }
  if (identical(interval, "Asymptotic")){
    output<-list(retMat,modelWeights,varMA)
    names(output)<-c("coef","weights","covar")
    return(invisible(output))
  }
  if (identical(interval, "Buckland")){
    output<-list(retMat,modelWeights)
    names(output)<-c("coef","weights")
    return(invisible(output))
  }
}
amavg<-function (modelList, modelWeights, expressions, interval=c("Buckland",
                                                                  "Asymptotic"),level=0.95, distribution="normal",seType="mod")
{
  require(sandwich, quietly = TRUE)
  require(car, quietly = TRUE)
  deltab<-function (object, g, func = g, ...)
  {
    if (!is.character(g))
      stop("The argument 'g' must be a character string")
    para <- coef(object)
    if(inherits(object,"lm")){
      para.names <- names(coef(object))
      para.names[1] <- gsub("\\(Intercept\\)", "Intercept", para.names[1])
      names(para) <- para.names
    }else{
      if(inherits(object,"drc")){
        coefVec<-coef(object)
        para.names<-sapply(strsplit(names(coefVec), ":"), "[[", 1)
        names(para)<-para.names
      }else{
        para.names <- names(para)
      }}
    g <- parse(text = g)
    q <- length(para)
    for (i in 1:q) {
      assign(names(para)[i], para[i])
    }
    gd <- rep(0,q)
    for (i in 1:q) {
      gd[i] <- eval(D(g, names(para)[i]))
    }
    gd
  }
  makeIIDdecomp <- function(modelObject,g)
  {
    numObsUsed <- ifelse(inherits(modelObject, "coxph"),
                         modelObject$n, ifelse(inherits(modelObject,"nls"),
                                               length(predict(modelObject)),ifelse(inherits(modelObject,"drc"),
                                                                                   length(predict(modelObject)), nrow(modelObject$model))))
    res.df<-df.residual(modelObject)
    db<-deltab(modelObject,g)
    iidVec0 <- db %*% bread(modelObject) %*% t(estfun(modelObject))
    moNAac <- modelObject$na.action
    numObs <- numObsUsed + length(moNAac)
    iidVec <- rep(0, numObs)
    if (!is.null(moNAac)) {
      iidVec[-moNAac] <- sqrt(numObs/numObsUsed) * iidVec0
    } else {
      iidVec <- iidVec0
    }
    list(iidVec = iidVec, numObsUsed = numObsUsed, numObs = numObs, res.df=res.df)
  }
  iidList <- list()
  numModels <- length(modelList)
  for(i in 1:numModels)
  {
    iidList[[i]]<- makeIIDdecomp(modelList[[i]], expressions[[i]])
  }
  iidresp <- matrix(as.vector(unlist(lapply(iidList, function(listElt) {
    listElt[[1]]}))), nrow = numModels, byrow = TRUE)
  numObsUsed <- as.vector(unlist(lapply(iidList, function(listElt) {
    listElt[[2]]})))
  res.df.vec <- as.vector(unlist(lapply(iidList, function(listElt) {
    listElt[[4]]})))
  thetaEst <- rep(NA, numModels)
  thetaSe <- rep(NA, numModels)
  for(i in 1:numModels)
  {
    if (inherits(modelList[[i]], "drc")){
      coefVec <- coef(modelList[[i]])
      names(coefVec) <- sapply(strsplit(names(coefVec), ":"), "[[", 1)
      deltaRes <- deltaMethod(coefVec,expressions[[i]],vcov(modelList[[i]]))
    } else {
      deltaRes <- deltaMethod(modelList[[i]],expressions[[i]])
    }
    thetaEst[i] <- deltaRes[1]
    thetaSe[i] <- deltaRes[2]
  }
  thetaEst <- unlist(thetaEst)
  thetaSe <- unlist(thetaSe)
  thetaMA <- as.numeric(modelWeights%*%thetaEst)
  if(identical(distribution, "t")){
    quant <- qt(1 - (1 - level)/2, df.residual(modelList[[1]]))
  } else {
    quant <- qnorm(1 - (1 - level)/2)
  }
  numObs <- iidList[[1]]$numObs
  covar <- (iidresp %*% t(iidresp)) / numObs
  vcMat <- covar / numObs # Defining the finite-sample variance-covariance matrix
  ## Replacing sandwich estimates by model-based standard errors
  modbas <- seType == "mod"
  if (any(modbas))
  {
    corMat <- cov2cor(vcMat)
    ## Retrieving standard errors for the specified estimate from the individual fits
    modSE <- thetaSe
    sanSE <- sqrt(diag(vcMat))
    sanSE[modbas] <- modSE[modbas]
    vcMat <- diag(sanSE,nrow=length(sanSE)) %*% corMat %*% diag(sanSE,nrow=length(sanSE))
  }
  if (identical(interval, "Asymptotic"))
  {
    numObs <- iidList[[1]]$numObs
    varMA <- modelWeights %*% t(modelWeights %*% vcMat)
    seMA <- sqrt(diag(varMA))
    quantVal <- quant * seMA
    retMat <- as.matrix(cbind(thetaMA, seMA, thetaMA - quantVal, thetaMA + quantVal))
    colnames(retMat) <- c("Estimate", "Std. Error", "Lower", "Upper")
  }
  if (identical(interval, "Buckland"))
  {
    seVec <- apply(sqrt(thetaSe^2 + (t(t(thetaEst) - mean(thetaMA)))^2) * modelWeights, 2,
                   sum)
    quantVal <- quant * seVec
    retMat <- as.matrix(cbind(thetaMA, seVec, thetaMA - quantVal, thetaMA + quantVal))
    colnames(retMat) <- c("Estimate", "Std. Error", "Lower", "Upper")
  }
  if (identical(interval, "Asymptotic")){
    output<-list(retMat,modelWeights,varMA)
    names(output)<-c("coef","weights","covar")
    return(invisible(output))
  }
  if (identical(interval, "Buckland")){
    output<-list(retMat,modelWeights)
    names(output)<-c("coef","weights")
    return(invisible(output))
  }
}
