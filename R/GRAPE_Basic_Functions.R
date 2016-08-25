



#' Make pairwise order representation of a sample
#'
#' Takes in a vector of gene expression values and returns a binary vector consisting of the pairwise rankings for the sample
#' @param samp A vector of gene expression values
#' @return Binary vector of the pairwise ranking representation of the samples
#' @examples
#' samp <- c(1,3,2,1.5)
#' makePairwiseOrder(samp)
#' @export
makePairwiseOrder <- function(samp){
  len <- length(samp)
  numpairs <- len*(len-1)/2
  pairwise_order <- vector(mode="numeric",numpairs)
  count <- 1
  for(j in 1:(len-1)){
    for(k in (j+1):len){
      if(samp[j] < samp[k]){pairwise_order[count] <- 1
      } else if(samp[j] > samp[k]) { pairwise_order[count] <- 0
      } else if(samp[j]==samp[k]){ ## decide tie arbitrary
        if(runif(1) > 0.5){pairwise_order[count] <- 1
        }else{pairwise_order[count] <- 0}
      }
      count <- count + 1
    }
  }
  return(pairwise_order)
}

#' Make template names from gene names
#'
#' Takes in vector of pathway gene names, returns names corresponding to the pairwise binary representation
#' @param path_genes A vector of pathway gene names
#' @return Names for the pairwise representation, of the form "gA<gB"
#' @examples
#' path_genes <- c("gene_A","gene_B","gene_C","gene_D")
#' makePairwiseOrderNames(path_genes)
#' @export
makePairwiseOrderNames <- function(path_genes){
  len <- length(path_genes)
  numpairs <- len*(len-1)/2
  pairnames <- rep(NA, numpairs)
  count <- 1
  for(j in 1:(len-1)){ for(k in (j+1):len){
    pairnames[count] <- paste0(path_genes[j]," < ",path_genes[k])
    count <- count + 1
  }}
  return(pairnames)
}

#' Make binary template and probability template
#'
#' Takes in matrix, where columns are samples and rows are pathway genes, outputs the binary and probability templates
#' @param submat A matrix where columns are samples and rows are pathway genes
#' @return List containing binary template vector and probability template vector
#' @examples
#' submat <- cbind(c(1,3,2,1.5),c(2,3,1.5,1.2),c(1.4,4.2,3.5,3.8))
#' rownames(submat) <- c("gene_A","gene_B","gene_C","gene_D")
#' temp <- makeBinaryTemplateAndProbabilityTemplate(submat)
#' bt <- temp$binary_template; pt <- temp$probability_template
#' cbind(bt,pt)
#' @export
makeBinaryTemplateAndProbabilityTemplate <- function(submat){
  ### Input: submat is matrix, columns are samples, rows are pathway genes
  ### Output: template is a list,
  ###         template$binary_template is binary template
  ###         template$probability_template is probability template
  if(is.null(rownames(submat))){print("Error: submat missing rownames.")}
  path_genes <- rownames(submat)
  len <- nrow(submat)
  numpairs <- len*(len-1)/2
  ### make pair names first
  pairnames <- makePairwiseOrderNames(path_genes)

  ### Make binary matrix of each sample
  binary_path_mat <- matrix(NA,nrow=numpairs,ncol=ncol(submat))
  rownames(binary_path_mat) <- pairnames
  colnames(binary_path_mat) <- colnames(submat)
  for(samp in 1:ncol(binary_path_mat)){
    binary_path_mat[,samp] <- makePairwiseOrder(submat[,samp])
  }
  probability_template <- rowMeans(binary_path_mat)
  binary_template <- rep(0,numpairs)
  names(binary_template) <- pairnames
  binary_template[which(probability_template > 0.5)] <- 1
  return(list(binary_template=binary_template,probability_template=probability_template))
}


