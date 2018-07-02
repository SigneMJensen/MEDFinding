POC<-function(object, data, type = "normal", pocMethods=c("linear","Dunnett","DRModels","MCP-MOD","user"), userContr=NA, useModels=NA){
  POC.models <- list()
  POC.tests <- list()
  
  if(type == "normal"){
    df <- data
    df[["dose"]] <- df[[ object[[3]] ]]
    df[["doseF"]] <- as.factor(df[["dose"]])
    df[["response"]] <- df[[ object[[2]] ]] 
    
    k<-0
    if("linear" %in% pocMethods){
      POC.models[[k+1]] <- lm(response ~ dose, data = df)  
      POC.tests[[k+1]] <- "dose"
      k<-k+1
    }
    if("DRModels" %in% pocMethods){
      j<-1
      for(i in 1:length(useModels)){
        options(show.error.messages = FALSE)
        POC.models[[k+j]] <- try(drm(response ~ dose, data = df, fct = useModels[[i]]),TRUE)
        if(!inherits(POC.models[[k+i]], "try-error")){
          POC.tests[[k+j]] <- "b"
          j <- j+1
        }
      }
      k <- k+j-1
    }
    
    if("Dunnett" %in% pocMethods){
      DModel <- lm(response ~ doseF -1, data = df)
      Dunnet.names <- paste("doseF",levels(df$doseF)[2:length(coef(m2))], 
                            " - doseF", levels(df$doseF)[1],sep="")
      
      for(i in 1:length(Dunnet.names)){
        POC.models[[k+i]] <- DModel
        POC.tests[[k+i]] <- Dunnet.names[i]
      }
    }
  }
  if(type=="binomial"){
    df <- data
    df[["dose"]] <- df[[ object[[3]] ]]
    df[["doseF"]] <- as.factor(df[["dose"]])
    df[["event"]] <- df[[ strsplit(as.character(object[[2]]), "/")[[2]] ]] 
    df[["total"]] <- df[[ strsplit(as.character(object[[2]]), "/")[[3]] ]] 
    
    if("linear" %in% pocMethods){
      k<-length(POC.tests)
      POC.models[[k+1]] <- glm(event/total ~ dose, data = df, weights = total, family = binomial)  
      POC.tests[[k+1]] <- "dose"
    }
    if("DRModels" %in% pocMethods){
      k<-length(POC.tests)
      j<-1
      for(i in 1:length(useModels)){
        #options(show.error.messages = FALSE)
        POC.models[[k+j]] <- try(drm(event/total ~ dose, data = df, type="binomial", weight=total, fct = useModels[[i]]),TRUE)
        if(!inherits(POC.models[[k+i]], "try-error")){
          POC.tests[[k+j]] <- "b"
          j <- j+1
        }
      }
    }
    if("Dunnett" %in% pocMethods){
      k<-length(POC.tests)
      DModel <- glm(event/total ~ doseF - 1, data = df, weights = total, family = binomial)
      Dunnet.names <- paste("doseF",levels(df$doseF)[2:length(coef(DModel))], 
                            " - doseF", levels(df$doseF)[1],sep="")
      
      for(i in 1:length(Dunnet.names)){
        POC.models[[k+i]] <- DModel
        POC.tests[[k+i]] <- Dunnet.names[i]
      }
    }
    if(("MCP-MOD" %in% pocMethods) | ("user" %in% pocMethods)){
      k<-length(POC.tests)
      contr<-list()
      AnovaModel <- glm(event/total ~ doseF - 1, data = df, weights = total, family = binomial)
      for(i in 1:ncol(userContr)){
        contr[[i]] <- paste(paste("(",userContr[,i],") * doseF",levels(df$doseF)[1:length(coef(DModel))],sep=""), collapse= " + ")
        POC.models[[k+i]] <- AnovaModel
        POC.tests[[k+i]] <- contr[[i]]
      }
    }
  }
  m <- mjust(POC.models,POC.tests, seType = "mod")  
  
  if("linear" %in% pocMethods){
    POC.pval <- summary(POC.models[[1]])$coef["dose","Pr(>|z|)"]
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
