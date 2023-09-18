#################################################################
# Filename: 	Wilks_importance.R
# Part of the SCE package, https://github.com/loong2020/Stepwise-Clustered-Ensemble.git
# Created: 		2019/05/17, Regina, SK, Canada
# Author: 		Kailong Li
# Email:		lkl98509509@gmail.com
# ===============================================================
Wilks_importance <- function(model,OOB_weight=FALSE)
{
  #: Extract information from Trees
  Wilk_mat <- lapply(model,function(x) data.frame(x$Tree))
  #: replace the wilks_min smaller than zero to zero
  Wilk_mat <- lapply(Wilk_mat,function(x) {x$wilk_min[which(x$wilk_min <= 0)] <- 0; return(x)})
  #: calculate the importance of each tree
  Wilk_mat <- lapply(Wilk_mat, function(x) data.frame(x,Imp=((x$lfMat+x$rtMat)/(x$lfMat[1]+x$rtMat[1]))*(1-x$wilk_min)))
  #: calculate the weighted contribution of each tree
  Imp <- lapply(Wilk_mat, function(x) aggregate(x, by=list(x$xCol), FUN = "sum"))
  #: give the true predictor index to importance of each tree
  Imp <- mapply(function(x,y) {x[-1,"Group.1"]<-y[x[-1,"Group.1"]];return(x[-1,])},x=Imp,y=lapply(model,function(x)x$Feature),SIMPLIFY=FALSE)
  if (OOB_weight==TRUE)
  {
    #: importance of each tree will be weighted by OOB error
    Imp_final <- mapply(function(x,y) data.frame(col_index=x[,"Group.1"],Importance=matrix(x[,"Imp"])%*%y),x=Imp,y=lapply(model,function(x)x$weight),SIMPLIFY=FALSE)
    #: Combine all trees to calculate the importance
    Imp_final <- do.call(rbind,Imp_final)
    Imp_final <- aggregate(Imp_final, by=list(Imp_final$col_index), FUN = "sum")[,-1]/table(do.call(c,lapply(model,function(x) x$Feature)))
    #: Extract the Feature name
    Feature_name <- unique(do.call(rbind,lapply(model,function(x) data.frame(name=x$XName,index=x$Feature))))
    Feature_name <- Feature_name[order(as.numeric(Feature_name$index)),"name"]
    #: Re-scale the Importance
    Imp_final[,-1] <- data.frame(apply(data.frame(Imp_final[,-1]),2,function(x)x/sum(x)))
    Imp_final$col_index <- Feature_name
    colnames(Imp_final) <- c("col_index",model[[1]]$YName)
    return(Imp_final)
  } else {
    Imp_final <- do.call(rbind,Imp)
    Imp_final <- aggregate(Imp_final, by=list(Imp_final$Group.1), FUN = "median")[,"Imp"]/table(do.call(c,lapply(model,function(x) x$Feature)))
    Imp_final <- data.frame(col_index=names(Imp_final),Importance=as.numeric(Imp_final))
    #: Extract the Feature name
    Feature_name <- unique(do.call(rbind,lapply(model,function(x) data.frame(name=x$XName,index=x$Feature))))
    Feature_name <- Feature_name[order(as.numeric(Feature_name$index)),"name"]
    #: Re-scale the Importance
    Imp_final[,-1] <- data.frame(apply(data.frame(Imp_final[,-1]),2,function(x)x/sum(x)))
    Imp_final$col_index <- Feature_name
    return(Imp_final)
  }
}