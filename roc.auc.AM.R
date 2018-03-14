roc.auc.AM <- function(df=NULL,fn=NULL,breaks,control,br.point=2,label,lGenes=FALSE){
  if(!is.null(fn)){
    x <- read.table(df,header = T, sep = "\t")
  } else {
    x = df
  }
  x$cut <- cut(x$FPKM,breaks,include.lowest = T, right = T)
  lCuts <- rev(levels(x$cut))
  mat <- matrix(nrow=7,ncol=length(lCuts),
                dimnames=list(c("normal","noise","sensitivity",
                                "1-specificity","specificity",
                                "Youden.Index","minDistance"),lCuts))
  if(lGenes){
    control<-tolower(control)
    idx <- which(tolower(x[,1]) %in% control)
    xtrue <- x[idx,]
    xnot <- x[-idx,]
    mat[1,] <- rev(table(xtrue$cut))
    mat[2,] <- rev(table(xnot$cut))
  } else {
      for( i in 1:br.point){
        genes1 <- unique(x[which(x$cut == lCuts[i]),1])
        genes2 <- unique(control[which(control$cut == lCuts[i]),1])
        mat[1,i] <- sum(genes1 %in% genes2)
        mat[2,i] <- sum(! genes1 %in% genes2)
      }
      for( i in (br.point+1):length(lCuts)){
        genes1 <- unique(x[which(x$cut == lCuts[i]),1])
        genes2 <- unique(control[which(control$cut == lCuts[i]),1])
        mat[2,i] <- sum(genes1 %in% genes2)
        mat[1,i] <- sum(! genes1 %in% genes2)
      }
  }
  TPR = c()
  FPR = c()
  spic = c()
  sumPos <- sum(mat[1,])
  sumNeg <- sum(mat[2,])
  for( i in seq(length(lCuts))){
    TPR <- c(TPR,round(sum(mat[1,1:i])/sumPos,2))
    FPR <- c(FPR,round(sum(mat[2,1:i])/sumNeg,2))
    spic <- c(spic,1-round(sum(mat[2,1:i])/sumNeg,2))
  }
  mat[3,] <- TPR
  mat[4,] <- FPR
  mat[5,] <- spic
  # calculate "Youden Index" J=max[Sn+Sp]
  mat[6,] <- TPR + spic
  # calculate distance from (0,1)
  mat[7,] <- round(sqrt((1-TPR)^2 + FPR^2),2)
  # The area under the ROC curve, when evaluated by the trapezoidal
  # rule, is given by
  mat.roc <- cbind(c(0,1),mat[c(3,5),])
  arr <- cbind(seq(ncol(mat.roc),2),seq((ncol(mat.roc)-1),1))
  auc <- sum(apply(arr,1,
                   FUN=function(x)
                     (mat.roc[1,x[1]]+mat.roc[1,x[2]])*
                     (mat.roc[2,x[2]]-mat.roc[2,x[1]])*.5))
  mat.roc[2,] <- 1-mat.roc[2,]
  rownames(mat.roc)[2] <- "1-specificity"
  return(list(auc=auc,label=label,br.point=lCuts[br.point],
              mat=mat,mat.roc=mat.roc))
}
