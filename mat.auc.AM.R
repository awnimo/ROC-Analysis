mat.auc.AM <- function(control=NULL,breaks,id=NULL,br.point=2,
                    depth,divide=100,genes=NULL,lGenes=FALSE){
  if(is.null(control)){
    ctrl <- genes
  } else {
    ctrl <- read.table(control, header = T, sep = "\t")
    # breaks
    ctrl$cut <- cut(ctrl$FPKM,breaks,include.lowest = T, right = T)
  }
  # observed data
  # filenames <- list.files(
  #   pattern = paste0(id,".*rep[1-4].genes.fpkm_tracking"))
  labels <- as.numeric(sub("x","",depth))/divide
  data <- list()
  for(i in 1:length(depth)){
    rep1 <- read.table(
      paste(id,depth[i],"rep1.genes.fpkm_tracking",sep="_"),
      header = T, sep = "\t", stringsAsFactors = F)[,c(7,1,10)]
    rep2 <- read.table(
      paste(id,depth[i],"rep2.genes.fpkm_tracking",sep="_"),
      header = T, sep = "\t", stringsAsFactors = F)[,c(7,1,10)]
    rep3 <- read.table(
      paste(id,depth[i],"rep3.genes.fpkm_tracking",sep="_"),
      header = T, sep = "\t", stringsAsFactors = F)[,c(7,1,10)]
    rep4 <- read.table(
      paste(id,depth[i],"rep4.genes.fpkm_tracking",sep="_"),
      header = T, sep = "\t", stringsAsFactors = F)[,c(7,1,10)]
    reps <- merge(rep1,rep2,by=c(1,2),suffixes = c(".rep1",".rep2"))
    reps <- merge(reps,rep3,by=c(1,2))
    reps <- merge(reps,rep4,by=c(1,2),suffixes = c(".rep3",".rep4"))
    reps$FPKM <- rowMeans(reps[,3:6])
    reps <- reps[,c(2,7)]
    d <- roc.auc(df=reps, breaks=breaks, control=ctrl,
                 br.point=br.point, label=labels[i], lGenes=lGenes)
    data[length(data)+1] <- list(d)
  }
  names(data) <- labels
  # plot auc vs depth
  matAUC <- matrix(nrow = length(data), ncol = 2, 
                   dimnames = list(c(),c("depth x 1e6","AUC")))
  for(i in 1:length(data)){
    matAUC[i,] <- c(data[[i]]$label,data[[i]]$auc)
  }
  plot(matAUC, type = "b")
  # plot ROC curves
  plot(data[[1]]$mat.roc[2,], data[[1]]$mat.roc[1,], type = "l", 
       col = 1, ylim = c(0,1), xlim = c(0,1), main = id,
       xlab = "1-specificity", ylab = "sensitivity")
  for(i in 2:length(data)){
    lines(data[[i]]$mat.roc[2,],data[[i]]$mat.roc[1,],type="l",col=i)
  }
  abline(a = 0, b = 1, lty = 3)
  return(list(name=id,mat.auc=matAUC,data=data))
}
