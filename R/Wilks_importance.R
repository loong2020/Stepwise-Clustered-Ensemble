#################################################################
# Filename: 	Wilks_importance.R
# Part of the SCE package, https://github.com/loong2020/Stepwise-Clustered-Ensemble.git
# Created: 		2019/05/17, Regina, SK, Canada
# Author: 		Kailong Li
# Email:		lkl98509509@gmail.com
# ===============================================================
Wilks_importance <- function(model,OOB_weight=TRUE)
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
  Imp <- mapply(function(x,y,z) {
    # Skip the root node (first row)
    if(nrow(x) <= 1) return(NULL)
    
    # Get the feature indices from the tree
    feature_indices <- x[-1, "Group.1"]
    
    # Create a mapping from feature indices to predictor names
    feature_map <- setNames(seq_along(y), y)
    
    # Map indices to predictor names
    predictor_names <- y[feature_indices]
    
    # Create new data frame with mapped indices
    result <- x[-1,]
    result$Group.1 <- predictor_names
    
    return(result)
  }, x=Imp, y=lapply(model,function(x)x$XName), z=lapply(model,function(x)x$Feature), SIMPLIFY=FALSE)
  
  # Remove NULL results
  Imp <- Imp[!sapply(Imp, is.null)]
  
  if (OOB_weight==TRUE)
  {
    #: importance of each tree will be weighted by OOB error
    Imp_final <- mapply(function(x,y) {
      # Ensure numeric values for importance calculation
      importance_values <- as.numeric(x[,"Imp"])
      weight <- as.numeric(y)
      data.frame(col_index=x[,"Group.1"], Importance=importance_values * weight)
    }, x=Imp, y=lapply(model,function(x)x$weight), SIMPLIFY=FALSE)
    
    #: Combine all trees to calculate the importance
    Imp_final <- do.call(rbind,Imp_final)
    Imp_final <- aggregate(Importance ~ col_index, data=Imp_final, FUN=sum)
    Imp_final$Importance <- Imp_final$Importance/length(model)
    
    #: Re-scale the Importance
    Imp_final$Importance <- Imp_final$Importance/sum(Imp_final$Importance)
    colnames(Imp_final) <- c("Predictor", "Relative_Importance")
    return(Imp_final)
  } else {
    Imp_final <- do.call(rbind,Imp)
    Imp_final <- aggregate(Importance ~ Group.1, data=data.frame(Group.1=Imp_final$Group.1, Importance=as.numeric(Imp_final$Imp)), FUN=median)
    Imp_final$Importance <- Imp_final$Importance/length(model)
    
    #: Re-scale the Importance
    Imp_final$Importance <- Imp_final$Importance/sum(Imp_final$Importance)
    colnames(Imp_final) <- c("Predictor", "Relative_Importance")
    return(Imp_final)
  }
}

#' Calculate importance scores for a single SCA tree model
#'
#' @param model A single SCA tree model
#' @return A data frame containing importance scores for each predictor
#' @export
SCA_importance <- function(model) {
  # Input validation
  if (!is.list(model)) {
    stop("model must be a SCA tree model")
  }
  
  # Extract tree information
  tree_df <- data.frame(model$Tree)
  
  # Replace wilks_min smaller than zero with zero
  tree_df$wilk_min[which(tree_df$wilk_min <= 0)] <- 0
  
  # Calculate importance scores
  tree_df$Imp <- ((tree_df$lfMat + tree_df$rtMat) / (tree_df$lfMat[1] + tree_df$rtMat[1])) * (1 - tree_df$wilk_min)
  
  # Aggregate importance by predictor
  imp_df <- aggregate(tree_df, by = list(tree_df$xCol), FUN = "sum")
  
  # Map feature indices to predictor names
  feature_indices <- imp_df[-1, "Group.1"]
  predictor_names <- model$XName[feature_indices]
  
  # Create result data frame
  result <- data.frame(
    Predictor = predictor_names,
    Raw_Importance = imp_df[-1, "Imp"],
    Split_Count = imp_df[-1, "lfMat"] + imp_df[-1, "rtMat"]
  )
  
  # Calculate relative importance
  result$Relative_Importance <- result$Raw_Importance / sum(result$Raw_Importance)
  
  return(result[,c("Predictor", "Relative_Importance")])
}