MEDFinding<-function(object, data, useModels, delta, type="normal",pocMethods,sMethod="AIC", userD = NA, 
                     margin=0.0001, userContr=NA, Adjusted=TRUE, Quantile=FALSE){
  fList <- mtList()
  if(type=="binomial"){
    df <- data
    df[["dose"]] <- df[[ object[[3]] ]]
    df[["doseF"]] <- as.factor(df[["dose"]])
    df[["event"]] <- df[[ strsplit(as.character(object[[2]]), "/")[[2]] ]] 
    df[["total"]] <- df[[ strsplit(as.character(object[[2]]), "/")[[3]] ]] 
    tmp <- POC(object, data=data, type = "binomial", pocMethods=c(pocMethods,"DRModels"), useModels=useModels,
               userContr=userContr)
    
    select.models <-tmp$Models
    select.tests <- tmp$Tests
    
    # Best model by AIC
    modelList<-list()
    for(i in 1:length(useModels)){
      #options(show.error.messages = FALSE)
      modelList[[i]] <- try(drm(event/total ~ dose, data = df, type="binomial", weight=total, fct = useModels[[i]]),TRUE)
      names(modelList)[[i]] <- paste("model.",i,sep='')
    }
    
    m.conv <- as.numeric(which(unlist(lapply(modelList, function(x) !inherits(x, "try-error")))))
    
    
    best.model <- names(which(unlist(lapply(modelList[m.conv],function(x) try(AIC(x),TRUE))) == 
                                min(unlist(lapply(modelList[m.conv],function(x) try(AIC(x),TRUE))))))
    
    if(sMethod=="user"){
      best.model <- names(modelList)[grepl(userD$name,unlist(lapply(modelList,function(x){x$fct$name})))]
    }
    slope <- ifelse(as.numeric(predict(modelList[[best.model]],data.frame(0))-predict(modelList[[best.model]],data.frame(Inf)))>0,
                    "decreasing",
                    "increasing")
    if( identical(slope,"increasing" )) {
      f0 <- ifelse(!is.na(coef(modelList[[best.model]])["c:(Intercept)"]), 
                   coef(modelList[[best.model]])["c:(Intercept)"], 
                   predict(modelList[[best.model]],data.frame(0)))
    } else {
      f0 <- ifelse(!is.na(coef(modelList[[best.model]])["d:(Intercept)"]), 
                   coef(modelList[[best.model]])["d:(Intercept)"], 
                   predict(modelList[[best.model]],data.frame(0)))
    } 
    crit1 <- ifelse(identical(slope,"increasing" ),f0+delta,f0-delta)
    
    # Test d-f0 = 0 for at tjekke om der i det hele taget findes in dosis der opfylder kravet
    # Find dosis0 hvor f(dosis0) = crit1 og derefter confidence interval for f(dosis0)
    if(sMethod %in% c("AIC","user")){
      dose.finding.models <- list()
      dose.finding.tests <- list()
      sim.adj <- list()
      model.exp <- fList[[modelList[[best.model]][["fct"]][["name"]]]]
      
      F.doses <- list() 
      F.doses[[1]] <- ED(modelList[[best.model]], crit1, type = "absolute", display = FALSE)[1]  
      
      lower <- -Inf
      upper <- Inf
      
      if( identical(slope,"increasing")) {
        dose.finding.models[[1]] <- modelList[[best.model]]
        dose.finding.tests[[1]] <- gsub("dose", paste(F.doses[[1]]) , paste(model.exp," - c",sep=""))
        
        m.dose.finding <- mjust(c(select.models,dose.finding.models),
                                c(select.tests,dose.finding.tests),
                                seType = "mod")
        
        beta <- m.dose.finding$coef[,"Estimate"]
        Sigma <- m.dose.finding$covar
        ifelse(Adjusted==TRUE,
               sim.adj[[1]]<- suppressWarnings(summary(glht(model = parm(beta, Sigma),alternative="greater"))),
               sim.adj[[1]]<- suppressWarnings(summary(glht(model = parm(beta, Sigma),alternative="greater"),test=adjusted(type="none"))))
        if(tail(sim.adj[[1]]$test$pvalues,1) < 0.05){ 
          lower <- F.doses[[1]]
          upper <- F.doses[[1]]
        } else {
          lower <- F.doses[[1]]
          upper <- Inf
        }
         if(upper-lower>margin){
        F.doses[[2]] <- max(df[["dose"]])
        dose.finding.models[[2]] <- modelList[[best.model]]
        dose.finding.tests[[2]] <- gsub("dose", paste(F.doses[[2]]) , paste(model.exp," - c",sep=""))
        
        m.dose.finding2 <- mjust(c(select.models,dose.finding.models),
                                 c(select.tests,dose.finding.tests),
                                 seType = "mod")
        
        beta2 <- m.dose.finding2$coef[,"Estimate"]
        Sigma2 <- m.dose.finding2$covar
        ifelse(Adjusted==TRUE,
               sim.adj[[2]]<- suppressWarnings(summary(glht(model = parm(beta2, Sigma2),alternative="greater"))),
               sim.adj[[2]]<- suppressWarnings(summary(glht(model = parm(beta2, Sigma2),alternative="greater"),test=adjusted(type="none"))))
        ifelse((tail(as.numeric(sim.adj[[2]]$test$coefficients),1)  >= delta) & 
                 (tail(sim.adj[[2]]$test$pvalues,1) < 0.05),  
               upper <- F.doses[[2]], upper <- 2)
         }
        adj.quantile<-list()
        
        l<-3
        while((upper-lower) > margin){
          F.doses[[l]] <- sum(c(lower,upper))/2
          dose.finding.models[[l]] <- modelList[[best.model]]
          dose.finding.tests[[l]] <- gsub("dose", paste(F.doses[[l]]) , paste(model.exp," - c",sep=""))
          
          m.dose.finding.a <- mjust(c(select.models,dose.finding.models),
                                    c(select.tests,dose.finding.tests),
                                    seType = "mod")
          beta.a <- m.dose.finding.a$coef[,"Estimate"]
          Sigma.a <- m.dose.finding.a$covar
          ifelse(Adjusted==TRUE,
                 sim.adj[[l]]<- suppressWarnings(summary(glht(model = parm(beta.a, Sigma.a),alternative="greater"))),
                 sim.adj[[l]]<- suppressWarnings(summary(glht(model = parm(beta.a, Sigma.a),alternative="greater"),test=adjusted(type="none"))))
          if(Quantile & Adjusted){
            adj.quantile[[l]] <- suppressWarnings(attr(confint(glht(model = parm(beta.a, Sigma.a),alternative="greater"))$confint,"calpha"))
          }
          if((tail(as.numeric(sim.adj[[l]]$test$coefficients),1)  >= delta) & (tail(sim.adj[[l]]$test$pvalues,1) < 0.05)){ 
            upper <- F.doses[[l]]
          } else {
            lower <- F.doses[[l]]
          } 
          l <- l+1  
        }
      }
    }
    if(sMethod=="MA"){
      myMod<-which(sapply(select.models,function(x)inherits(x,"drc")))
      as.numeric("linear" %in% pocMethods)
      maModelList <- list()
      maTestMat <- matrix(0,ncol = as.numeric("linear" %in% pocMethods) + 
                            length(useModels)+
                            as.numeric(sum(as.numeric(c("Dunnett","MCP-MOD","user")  %in% pocMethods))>0),
                          nrow = length(select.models))
      w.mat <- maTestMat
      k<-0
      if("linear" %in% pocMethods){
        maModelList[[k+1]] <- select.models[[1]] 
        maTestMat[,k+1] <- select.tests[[1]]
        w.mat[1,k+1] <- 1
      }
      k<-length(maModelList)
      for(i in 1:length(myMod)){
        maModelList[[k+i]] <- select.models[[myMod[i]]]
        maTestMat[,k+i] <- select.tests[[myMod[i]]]
        w.mat[k+i,k+i] <- 1
      }
      k<-length(maModelList)
      maModelList[[k+1]] <- select.models[[k+1]]
      if((k)<length(select.tests)){
        maTestMat[,k+1] <- select.tests[[k+1]]
        for(i in (k+1):length(select.tests)){
          maTestMat[i,k+1] <- select.tests[[i]]
          w.mat[i,k+1] <- 1
        }
      }
      
      w.AIC<-sapply(select.models[myMod],function(x){AIC(x)})
      sim.ma.adj<-list()
      
      model.exp <- lapply(useModels,function(x){fList[[x$name]]})
      
      F.doses <- list()
      F.doses[[1]]<-0.001
      ma1<- maTestMat[1,]
      ma1[myMod] <- sapply(model.exp,function(x){gsub("dose", paste(F.doses[[1]]) , paste(x," - c",sep=""))})
      maTestMat<-rbind(maTestMat,ma1)
      w1<-rep(0,ncol(w.mat))
      w1[myMod]<- exp((min(w.AIC)-w.AIC)/2)/sum(exp((min(w.AIC)-w.AIC)/2))
      w.mat<-rbind(w.mat,w1)
      
      weightMat <- matrix(0,ncol=length(maModelList)*nrow(maTestMat),nrow=nrow(maTestMat))
      for(i in 1:nrow(maTestMat)){
        weightMat[i,(i-1)*length(maModelList)+1:length(maModelList)]<-w.mat[i,]
      }
      ma<- amavg(rep(maModelList,nrow(maTestMat)), weightMat,as.list(t(maTestMat)),interval="Asymptotic")
      
      ifelse(Adjusted,
             sim.ma.adj[[1]] <- suppressWarnings(summary(glht(model = parm(ma$coef[,1], ma$covar),alternative="greater"))),
             sim.ma.adj[[1]] <- suppressWarnings(summary(glht(model = parm(ma$coef[,1], ma$covar),alternative="greater"),test=adjusted(type="none"))))
      
      if((tail(as.numeric(sim.ma.adj[[1]]$test$coefficients),1)  >= delta) & (tail(sim.ma.adj[[1]]$test$pvalues,1) < 0.05)){ 
        stop("Minimum effective dose < 0.001")
      } else {
        lower <- F.doses[[1]]
      } 
      
      F.doses[[2]] <- max(df[["dose"]])
      ma1[myMod] <- sapply(model.exp,function(x){gsub("dose", paste(F.doses[[2]]) , paste(x," - c",sep=""))})
      maTestMat<-rbind(maTestMat,ma1)
      w.mat<-rbind(w.mat,w1)
      
      weightMat <- matrix(0,ncol=length(maModelList)*nrow(maTestMat),nrow=nrow(maTestMat))
      for(i in 1:nrow(maTestMat)){
        weightMat[i,(i-1)*length(maModelList)+1:length(maModelList)]<-w.mat[i,]
      }
      ma<- amavg(rep(maModelList,nrow(maTestMat)), weightMat,as.list(t(maTestMat)),interval="Asymptotic")
      
      ifelse(Adjusted,
             sim.ma.adj[[2]] <- suppressWarnings(summary(glht(model = parm(ma$coef[,1], ma$covar),alternative="greater"))),
             sim.ma.adj[[2]] <- suppressWarnings(summary(glht(model = parm(ma$coef[,1], ma$covar),alternative="greater"),test=adjusted(type="none"))))
      
      if((tail(as.numeric(sim.ma.adj[[2]]$test$coefficients),1)  >= delta) & (tail(sim.ma.adj[[2]]$test$pvalues,1) < 0.05)){ 
        upper <- F.doses[[2]]
      } else {
        stop("Minimum effective dose > largest dose")
      } 
      
      adj.quantile<-list()
      
      l<-3
      while((upper-lower) > margin){
        F.doses[[l]] <- sum(c(lower,upper))/2
        ma1[myMod] <- sapply(model.exp,function(x){gsub("dose", paste(F.doses[[l]]) , paste(x," - c",sep=""))})
        maTestMat<-rbind(maTestMat,ma1)
        w.mat<-rbind(w.mat,w1)
        
        weightMat <- matrix(0,ncol=length(maModelList)*nrow(maTestMat),nrow=nrow(maTestMat))
        for(i in 1:nrow(maTestMat)){
          weightMat[i,(i-1)*length(maModelList)+1:length(maModelList)]<-w.mat[i,]
        }
        ma<- amavg(rep(maModelList,nrow(maTestMat)), weightMat,as.list(t(maTestMat)),interval="Asymptotic")
        
        ifelse(Adjusted,
               sim.ma.adj[[l]] <- suppressWarnings(summary(glht(model = parm(ma$coef[,1], ma$covar),alternative="greater"))),
               sim.ma.adj[[l]] <- suppressWarnings(summary(glht(model = parm(ma$coef[,1], ma$covar),alternative="greater"),test=adjusted(type="none"))))
        if(Quantile & Adjusted){
          adj.quantile[[l]] <- suppressWarnings(attr(confint(glht(model = parm(ma$coef[,1], ma$covar),alternative="greater"))$confint,"calpha"))
        }
        if((tail(as.numeric(sim.ma.adj[[l]]$test$coefficients),1)  >= delta) & (tail(sim.ma.adj[[l]]$test$pvalues,1) < 0.05)){ 
          upper <- F.doses[[l]]
        } else {
          lower <- F.doses[[l]]
        } 
        l <- l+1  
      }
    }
  }
  
  MEDose<-list(MED = tail(F.doses,1),
               AllDoses = F.doses,
               size = length(F.doses), 
               Quantile = tail(adj.quantile,1),
               All.quantiles = adj.quantile,
               SIM = sim.ma.adj)
  class(MEDose) = "MEDFinding"
  MEDose
}
