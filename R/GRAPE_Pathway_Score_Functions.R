######################################################################
### Pathway Score Functions
######################################################################

#' Calculate Pathway Scores
#'
#' Calculate pathway scores of a single pathway of a set of samples relative to a reference set of samples
#'
#' @param refmat Pathway expression matrix of reference samples. Rows are genes, columns are samples.
#' @param newmat Pathway expression matrix of new samples. Rows are genes, columns are samples.
#' @param w Weight function. Default is quadratic weight function.
#' @return Vector of pathway scores of each sample in newmat.
#' @importFrom stats median
#' @importFrom stats quantile
#' @importFrom stats runif
#' @examples
#' ## Toy example: 50 reference samples
#' set.seed(10);
#' refmat <- matrix(rnorm(5*50),nrow=5,ncol=50); rownames(refmat) <- paste0("g",1:5)
#' ### make g2 and g5 larger in refmat
#' refmat[2,] <- rnorm(50,3,2); refmat[5,] <- rnorm(50,4,4)
#' ### 15 new samples
#' newmat <- matrix(rnorm(5*15),nrow=5,ncol=15); rownames(newmat) <- paste0("g",1:5)
#' ### make g2 and g3 larger in newmat
#' newmat[2,] <- rnorm(15,2,3); newmat[3,] <- rnorm(15,4,3)
#' ps_new <- getPathwayScores(refmat,newmat) ### get pathway scores of new samples
#' ps_ref <- getPathwayScores(refmat,refmat) ### get pathway scores of reference samples
#' ps_both <- getPathwayScores(refmat,cbind(refmat,newmat)) ### get pathway scores of both
#' # > ps_new
#' # [1]  6.2720  8.5696  9.9904  6.9056  3.7824  8.9344 13.0880 10.2912  3.7824
#' # 0.0384 13.1136  6.8032  4.8512 12.8512 10.2912
#' @export
getPathwayScores <- function(refmat,newmat,w=w_quad){
  ### Input: refmat is matrix of n reference samples (columns) and m pathway genes
  ### newmat is matrix of n non-reference (new) samples (columns) and m pathway genes
  ### w is weight function
  ### Output: pathway_scores is vector of the pathway scores of each new sample relative to the reference samples
  temp <- makeBinaryTemplateAndProbabilityTemplate(refmat)
  bt <- temp$binary_template
  pt <- temp$probability_template
  denominator <- sum(w(pt))

  ### Get distances of reference samples
  ref_dists <- rep(NA,ncol(refmat))
  for(smp in 1:ncol(refmat)){
    pwo <- makePairwiseOrder(refmat[,smp])
    missidc <- which(abs(pwo-bt)==1)
    goodidc <- setdiff(1:length(pwo),missidc)
    if(length(goodidc)==0)ref_dists[smp] <- 0 ### avoid error in case of complete mismatch
    if(length(goodidc>0)){
      accuracy <- sum(w(pt[goodidc]))/denominator
      ref_dists[smp] <- 100*(1-accuracy)
    }
  }
  ### Get distances of non-reference samples
  new_dists <- rep(NA,ncol(newmat))
  for(smp in 1:ncol(newmat)){
    pwo <- makePairwiseOrder(newmat[,smp])
    missidc <- which(abs(pwo-bt)==1)
    goodidc <- setdiff(1:length(pwo),missidc)
    if(length(goodidc)==0)new_dists[smp] <- 0 ### avoid error in case of complete mismatch
    if(length(goodidc>0)){
      accuracy <- sum(w(pt[goodidc]))/denominator
      new_dists[smp] <- 100*(1-accuracy)
    }
  }

  theta <- median(ref_dists) ### median
  q <- quantile(ref_dists)
  delta <- q[4] - q[2] ### inter-quartile distance

  ### modification in case delta is too small to avoid infinity error
  denominator <- sum(w(pt))
  ### min iqd is set to be a single violation
  miniqd <- 100*w(.75)/denominator

  delta <- max(delta,miniqd)

  pathway_scores <- (new_dists-theta)/delta
  pathway_scores[which(pathway_scores < 0)] <- 0
  return(pathway_scores)
}

#' Calculate Pathway Space Matrix
#'
#' Represents new samples as vectors of pathway scores relative to reference samples
#'
#' @param refge Gene expression matrix of reference samples. Rows are genes, columns are samples.
#' @param newge Gene expression matrix of new samples. Rows are genes, columns are samples.
#' @param pathway_list List of pathways. Each pathway is a character vector consisting of gene names.
#' @param w Weight function. Default is quadratic weight function.
#' @return Vector of pathway scores of each sample in newmat.
#' @examples
#' #' ### Make pathway scores mat
#' set.seed(10)
#' ### 50 reference samples
#' refge <- matrix(rnorm(10*50),nrow=10,ncol=50); rownames(refge) <- paste0("g",1:10)
#' refge[c(2,5,8),] <- matrix(rnorm(3*50,mean=2,sd=2))
#' refge[c(3,4,7),] <- matrix(rnorm(3*50,mean=4,sd=4))
#' ### 6 new samples
#' newge <- matrix(rnorm(10*6),nrow=10,ncol=6); rownames(newge) <- paste0("g",1:10)
#' newge[c(2:7),] <- matrix(rnorm(6*6,mean=3,sd=1))
#' newge[c(1,9),] <- matrix(rnorm(2*6,mean=5,sd=3))
#' pathway_list <- list(set1=paste0("g",1:4),set2=paste0("g",5:10),set3=paste0("g",c(1,4,8:10)))
#' psmat <- makeGRAPE_psMat(refge,newge,pathway_list)
#' # > psmat
#' # [,1]     [,2]     [,3]     [,4]     [,5]     [,6]
#' # set1 2.397426 1.406275 2.516492 2.358809 2.555109 2.358809
#' # set2 0.670354 3.245575 3.962389 2.670354 1.741150 1.579646
#' # set3 1.536017 2.167373 2.167373 2.167373 2.148305 1.809322
#' @export
makeGRAPE_psMat <- function(refge,newge,pathway_list,w=w_quad){
  if(!identical(rownames(refge),rownames(newge))){print("Warning: different row names for refge and newge.")}
  npaths <- length(pathway_list)
  psmat <- matrix(NA,nrow=npaths,ncol=ncol(newge))
  rownames(psmat) <- names(pathway_list)
  for(pth in 1:npaths){
    path_genes <- pathway_list[[pth]]
    refmat <- refge[which(rownames(refge)%in%path_genes),]
    newmat <- newge[which(rownames(newge)%in%path_genes),]
    psmat[pth,] <- getPathwayScores(refmat,newmat,w)
  }
  return(psmat)
}
