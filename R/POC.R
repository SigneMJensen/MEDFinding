POC<-function(object, data, type = "normal", pocMethods=c("linear","Dunnett","DRModels","MCP-MOD","user"), userContr=NA, useModels=NA){
  fList<-mtList()
  POC.models <- list()
  POC.tests <- list()
  
  if(type == "normal"){
    df1 <- data
    df1[["dose"]] <- df1[[ object[[3]] ]]
    df1[["doseF"]] <- as.factor(df1[["dose"]])
    df1[["response"]] <- df1[[ object[[2]] ]] 
    
    k<-0
    if("linear" %in% pocMethods){
      POC.models[[k+1]] <- lm(response ~ dose, data = df1)  
      POC.tests[[k+1]] <- "dose"
      k<-k+1
    }
    if("DRModels" %in% pocMethods){
      j<-1
      for(i in 1:length(useModels)){
        options(show.error.messages = FALSE)
        POC.models[[k+j]] <- try(drm(response ~ dose, data = df1, fct = useModels[[i]]),TRUE)
        if(!inherits(POC.models[[k+i]], "try-error")){
          POC.tests[[k+j]] <- fList[[useModels[[i]]$name]]$poc
          j <- j+1
        }
      }
      k <- k+j-1
    }
    
    if("Dunnett" %in% pocMethods){
      DModel <- lm(response ~ doseF -1, data = df1)
      Dunnet.names <- paste("doseF",levels(df1$doseF)[2:length(coef(DModel))], 
                            " - doseF", levels(df1$doseF)[1],sep="")
      
      for(i in 1:length(Dunnet.names)){
        POC.models[[k+i]] <- DModel
        POC.tests[[k+i]] <- Dunnet.names[i]
      }
    }
  }
  if(type=="binomial"){
    df1 <- data
    df1[["dose"]] <- df1[[ object[[3]] ]]
    df1[["doseF"]] <- as.factor(df1[["dose"]])
    df1[["event"]] <- df1[[ strsplit(as.character(object[[2]]), "/")[[2]] ]] 
    df1[["total"]] <- df1[[ strsplit(as.character(object[[2]]), "/")[[3]] ]] 
    
    if("linear" %in% pocMethods){
      k<-length(POC.tests)
      POC.models[[k+1]] <- glm(event/total ~ dose, data = df1, weights = total, family = binomial)  
      POC.tests[[k+1]] <- "dose"
    }
    if("DRModels" %in% pocMethods){
      k<-length(POC.tests)
      j<-1
      for(i in 1:length(useModels)){
        #options(show.error.messages = FALSE)
        POC.models[[k+j]] <- try(drm(event/total ~ dose, data = df1, type="binomial", weight=total, fct = useModels[[i]]),TRUE)
        if(!inherits(POC.models[[k+i]], "try-error")){
          POC.tests[[k+j]] <- "b"
          j <- j+1
        }
      }
    }
    if("Dunnett" %in% pocMethods){
      k<-length(POC.tests)
      DModel <- glm(event/total ~ doseF - 1, data = df1, weights = total, family = binomial)
      Dunnet.names <- paste("doseF",levels(df1$doseF)[2:length(coef(DModel))], 
                            " - doseF", levels(df1$doseF)[1],sep="")
      
      for(i in 1:length(Dunnet.names)){
        POC.models[[k+i]] <- DModel
        POC.tests[[k+i]] <- Dunnet.names[i]
      }
    }
    if(("MCP-MOD" %in% pocMethods) | ("user" %in% pocMethods)){
      k<-length(POC.tests)
      contr<-list()
      AnovaModel <- glm(event/total ~ doseF - 1, data = df1, weights = total, family = binomial)
      for(i in 1:ncol(userContr)){
        contr[[i]] <- paste(paste("(",userContr[,i],") * doseF",levels(df1$doseF)[1:length(coef(DModel))],sep=""), collapse= " + ")
        POC.models[[k+i]] <- AnovaModel
        POC.tests[[k+i]] <- contr[[i]]
      }
    }
  }
  m <- mjust(POC.models,POC.tests, seType = "mod")  
  
  if("linear" %in% pocMethods){
    POC.pval <- tail(summary(POC.models[[1]])$coef["dose",],1)
  }
  if(length(POC.tests)>1){
    POC.pval<-suppressWarnings(summary(glht(model=parm(m$coef[,"Estimate"],m$covar),alternative="greater"))$test$pvalues)
  }
  results <- list(POC.TF = paste("Proof of concept: ", sum(POC.pval<0.05) > 0,sep=""),
                  Models = POC.models,
                  Tests = POC.tests,
                  Pval = POC.pval)
  #POC.pval
  class(results) = "POC"
  results
}
