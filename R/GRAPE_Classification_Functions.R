######################################################################
### GRAPE Classification
######################################################################
### GRAPE quadratic weight function
### x can be scalar or vecor

#' Quadratic weight function
#'
#' Calculates the weights of all input entries. All entries should take values in [0,1].
#' @param x Any number, vector of matrix.
#' @return Weight of each element
#' @examples
#' w_quad(0.95)
#' w_quad(cbind(c(.7,.8),c(.9,.1)))
#' @export
w_quad <- function(x){return(4*abs(x-0.5)^2)}

#' GRAPE Classification
#'
#' Classification of a samples according to grape distances from templates. Usually applied to the gene expression values for a single pathway.
#' @param trainmat Matrix of gene expression for set of genes accross training set samples. Each column is a sample.
#' @param testmat Matrix of gene expression for set of genes accross test set samples. Each column is a sample.
#' @param train_labels Vector of class labels for each sample in the training set.
#' @param w Weight function. Default is quadratic weight function.
#' @return Predicted class labels for test set
#' @examples
#' # Toy example of two classes
#' set.seed(10); path_genes <- c("gA","gB","gC","gD"); nsamps <- 50 # Four genes, 50 samples per class
#' class_one_samps <- matrix(NA,nrow=length(path_genes),ncol=nsamps) # Class 1
#' rownames(class_one_samps) <- path_genes
#' class_one_samps[1,] <- rnorm(ncol(class_one_samps),4,2)
#' class_one_samps[2,] <- rnorm(ncol(class_one_samps),5,4)
#' class_one_samps[3,] <- rnorm(ncol(class_one_samps),1,1)
#' class_one_samps[4,] <- rnorm(ncol(class_one_samps),2,1)
#' class_two_samps <- matrix(NA,nrow=length(path_genes),ncol=nsamps) # Class 2
#' rownames(class_two_samps) <- path_genes
#' class_two_samps[1,] <- rnorm(ncol(class_two_samps),2,3)
#' class_two_samps[2,] <- rnorm(ncol(class_two_samps),5,2)
#' class_two_samps[3,] <- rnorm(ncol(class_two_samps),1,1)
#' class_two_samps[4,] <- rnorm(ncol(class_two_samps),0,1)
#' all_samps <- cbind(class_one_samps,class_two_samps)
#' labels <- c(rep(1,nsamps),rep(2,nsamps))
#' testid <- sample.int(100,20)
#' trainmat <- all_samps[,-testid]
#' train_labels <- labels[-testid]
#' testmat <- all_samps[,testid]
#' test_labels <- labels[testid]
#' yhat <- predictClassGRAPE(trainmat,testmat,train_labels,w_quad)
#' sum(diag(table(test_labels,yhat)))/length(test_labels) # accuracy
#' # [1] 0.8
#' @export
predictClassGRAPE <- function(trainmat, testmat, train_labels, w=w_quad){
  ### Input:
  ### trainmat and testmat are matrices, columns are samples, rows are pathway genes
  ### train_labels is a vector of class labels for the training set
  ### w is weight function
  ### Output: yhat is vector of class membership predicted by GRAPE
  if(!identical(rownames(trainmat),rownames(testmat))){print("Error: non-identical features")}
  ### Make Templates and PT using train data
  num_genes <- nrow(trainmat)
  num_pairs <- num_genes*(num_genes-1)/2
  num_classes <- length(unique(train_labels))
  classes <- sort(unique(train_labels))

  bin_templates <- matrix(NA,nrow=num_pairs,ncol=num_classes)
  prob_templates <- matrix(NA,nrow=num_pairs,ncol=num_classes)
  colnames(bin_templates) <- classes
  colnames(prob_templates) <- classes
  for(j in 1:num_classes){
    submat <- trainmat[,which(train_labels==classes[j])]
    temp <- makeBinaryTemplateAndProbabilityTemplate(submat)
    bin_templates[,j] <- temp$binary_template
    prob_templates[,j] <- temp$probability_template
  }
  ### we have made binary templates and probability templates for all classes
  yhat <- rep(NA, ncol(testmat)) # Vector of predicted class labels

  ### go through each sample in the test set
  for(smp in 1:ncol(testmat)){
    pwo <- makePairwiseOrder(testmat[,smp])
    scores <- rep(NA,num_classes)
    for(class in 1:num_classes){
      bt <- bin_templates[,class]
      pt <- prob_templates[,class]
      denominator <- sum(w(pt))
      missidc <- which(abs(pwo-bt)==1)
      scores[class] <- sum(w(pt[missidc]))/denominator
    }
    yhat[smp] <- classes[which.min(scores)]
  }
  return(yhat)
}
######################################################################
### DIRAC Classification
######################################################################
#' DIRAC Classification
#'
#' Classification of a samples according to dirac distances from templates. Usually applied to the gene expression values for a single pathway.
#' @param trainmat Matrix of gene expression for set of genes accross training set samples. Each column is a sample.
#' @param testmat Matrix of gene expression for set of genes accross test set samples. Each column is a sample.
#' @param train_labels Vector of class labels for each sample in the training set.
#' @return Predicted class labels for test set
#' @examples
#' # Toy example of two classes
#' set.seed(10); path_genes <- c("gA","gB","gC","gD"); nsamps <- 50 # Four genes, 50 samples per class
#' class_one_samps <- matrix(NA,nrow=length(path_genes),ncol=nsamps) # Class 1
#' rownames(class_one_samps) <- path_genes
#' class_one_samps[1,] <- rnorm(ncol(class_one_samps),4,2)
#' class_one_samps[2,] <- rnorm(ncol(class_one_samps),5,4)
#' class_one_samps[3,] <- rnorm(ncol(class_one_samps),1,1)
#' class_one_samps[4,] <- rnorm(ncol(class_one_samps),2,1)
#' class_two_samps <- matrix(NA,nrow=length(path_genes),ncol=nsamps) # Class 2
#' rownames(class_two_samps) <- path_genes
#' class_two_samps[1,] <- rnorm(ncol(class_two_samps),2,3)
#' class_two_samps[2,] <- rnorm(ncol(class_two_samps),5,2)
#' class_two_samps[3,] <- rnorm(ncol(class_two_samps),1,1)
#' class_two_samps[4,] <- rnorm(ncol(class_two_samps),0,1)
#' all_samps <- cbind(class_one_samps,class_two_samps)
#' labels <- c(rep(1,nsamps),rep(2,nsamps))
#' testid <- sample.int(100,20)
#' trainmat <- all_samps[,-testid]
#' train_labels <- labels[-testid]
#' testmat <- all_samps[,testid]
#' test_labels <- labels[testid]
#' yhat <- predictClassDIRAC(trainmat,testmat,train_labels)
#' sum(diag(table(test_labels,yhat)))/length(test_labels) # accuracy
#' # [1] 0.7
#' @export
predictClassDIRAC <- function(trainmat, testmat, train_labels){
  ### Input:
  ### trainmat and testmat are matrices, columns are samples, rows are pathway genes
  ### train_labels is a vector of class labels for the training set
  ### Output: yhat is vector of class membership predicted by DIRAC
  if(!identical(rownames(trainmat),rownames(testmat))){print("Error: non-identical features")}
  ### Make Templates and PT using train data
  num_genes <- nrow(trainmat)
  num_pairs <- num_genes*(num_genes-1)/2
  num_classes <- length(unique(train_labels))
  classes <- sort(unique(train_labels))

  bin_templates <- matrix(NA,nrow=num_pairs,ncol=num_classes)
  colnames(bin_templates) <- classes
  for(j in 1:num_classes){
    submat <- trainmat[,which(train_labels==classes[j])]
    temp <- makeBinaryTemplateAndProbabilityTemplate(submat)
    bin_templates[,j] <- temp$binary_template
  }
  ### we have made binary templates for all classes
  yhat <- rep(NA, ncol(testmat)) # Vector of predicted class labels

  ### go through each sample in the test set
  for(smp in 1:ncol(testmat)){
    pwo <- makePairwiseOrder(testmat[,smp])
    scores <- rep(NA,num_classes)
    for(class in 1:num_classes){
      bt <- bin_templates[,class]
      scores[class] <- sum(abs(pwo - bt))/length(bt)
    }
    yhat[smp] <- classes[which.min(scores)]
  }
  return(yhat)
}
######################################################################
### Pathway Centroid (PC) Classification
######################################################################
#' PC Classification
#'
#' Classification of a samples according to euclidean distances from PC templates. Usually applied to the gene expression values for a single pathway.
#' @param trainmat Matrix of gene expression for set of genes accross training set samples. Each column is a sample.
#' @param testmat Matrix of gene expression for set of genes accross test set samples. Each column is a sample.
#' @param train_labels Vector of class labels for each sample in the training set.
#' @return Predicted class labels for test set
#' @examples
#' # Toy example of two classes
#' set.seed(10); path_genes <- c("gA","gB","gC","gD"); nsamps <- 50 # Four genes, 50 samples per class
#' class_one_samps <- matrix(NA,nrow=length(path_genes),ncol=nsamps) # Class 1
#' rownames(class_one_samps) <- path_genes
#' class_one_samps[1,] <- rnorm(ncol(class_one_samps),4,2)
#' class_one_samps[2,] <- rnorm(ncol(class_one_samps),5,4)
#' class_one_samps[3,] <- rnorm(ncol(class_one_samps),1,1)
#' class_one_samps[4,] <- rnorm(ncol(class_one_samps),2,1)
#' class_two_samps <- matrix(NA,nrow=length(path_genes),ncol=nsamps) # Class 2
#' rownames(class_two_samps) <- path_genes
#' class_two_samps[1,] <- rnorm(ncol(class_two_samps),2,3)
#' class_two_samps[2,] <- rnorm(ncol(class_two_samps),5,2)
#' class_two_samps[3,] <- rnorm(ncol(class_two_samps),1,1)
#' class_two_samps[4,] <- rnorm(ncol(class_two_samps),0,1)
#' all_samps <- cbind(class_one_samps,class_two_samps)
#' labels <- c(rep(1,nsamps),rep(2,nsamps))
#' testid <- sample.int(100,20)
#' trainmat <- all_samps[,-testid]
#' train_labels <- labels[-testid]
#' testmat <- all_samps[,testid]
#' test_labels <- labels[testid]
#' yhat <- predictClassPC(trainmat,testmat,train_labels)
#' sum(diag(table(test_labels,yhat)))/length(test_labels) # accuracy
#' # [1] 0.55
#' @export
predictClassPC <- function(trainmat, testmat, train_labels){
  ### Input:
  ### trainmat and testmat are matrices, columns are samples, rows are pathway genes
  ### train_labels is a vector of class labels for the training set
  ### Output: yhat is vector of class membership predicted by PC
  if(!identical(rownames(trainmat),rownames(testmat))){print("Error: non-identical features")}
  ### Make Templates and PT using train data
  num_genes <- nrow(trainmat)
  num_classes <- length(unique(train_labels))
  classes <- sort(unique(train_labels))

  ### Make centroids
  centroids <- matrix(NA,nrow=nrow(trainmat),ncol=num_classes)
  colnames(centroids) <- classes
  for(j in 1:num_classes){
    submat <- trainmat[,which(train_labels==classes[j])]
    centroids[,j] <- rowMeans(submat)
  }
  ### we have centroids for all classes
  yhat <- rep(NA, ncol(testmat)) # Vector of predicted class labels

  ### go through each sample in the
  for(smp in 1:ncol(testmat)){
    samp <- testmat[,smp]
    ### dirac
    scores <- rep(NA,num_classes)
    for(class in 1:num_classes){scores[class] <- sqrt(sum((centroids[,class]-samp)^2))}
    yhat[smp] <- classes[which.min(scores)]
  }
  return(yhat)
}
