#' @title Edgewise Linear Model
#'
#' @description This function computes the specified linear model (LM) for each edge in the network,
#'  and calculates the multiple testing p-value based on the p.adjust function.
#'
#' @usage edge_lm(net, nnodes, idata, mod, diag = FALSE, padj,
#'        cores = NULL, expList = NULL, verbose = TRUE,
#'        ...)
#'
#' @param net 3D volume (2D matrices for each observation) or 2D matrix of edges as columns.
#' @param nnodes Number of network nodes.
#' @param idata Matrix or data.frame including independent variables of interest of the model.
#' @param mod Model, specify as a string, e.g., "~Group + Age".
#' @param diag Logical indicating if matrix diagonal is to be included in the analysis (default: FALSE).
#' @param padj Character string that indicates the p.adjust method.
#' @param cores Number of selected cores for parallel computing (default: NULL).
#' @param expList Character string adding variable names to the varlist of 'clusterExport' (default: NULL).
#' @param verbose Logical indicating if messages should be printed (default: TRUE).
#' @param ... Additional arguments to be passed to the low level 'lm' function.
#'
#' @details It's VERY IMPORTANT when giving \emph{net} as a 2D matrix or data.frame, to be
#'  completely sure that column distribution fits that of the upper triangle indices of an
#'  \emph{nnodes} * \emph{nnodes} matrix. This may be verified through the edge indices, e.g.,
#'  "which(upper.tri(matrix(nrow = nnodes, ncol = nnodes)), arr.ind = T)" (see vignette NBR-LME
#'  for more details).
#'
#' To know more about \emph{padj} methods, check help for the \code{\link{p.adjust}} function. It is
#' noticeable that this multiple comparison approach can be much more faster than the permutations run by
#' the Network-Based Statistics framework, however this is a much more conservative approach
#' (see Zalesky et al. (2010) \url{https://doi.org/10.1016/j.neuroimage.2010.06.041} for more detail).
#'
#' @return \emph{data.frame} containing the edge labels, observed statistics,
#'  their corresponding p-value, and their adjusted p-values.
#'
#' @examples
#' data(frontal2D)
#' \donttest{
#' edge_result <- edge_lm(net = frontal2D[,-(1:3)], nnodes = 28,
#'   idata = frontal2D[,1:3], mod = "~ Group + Sex * Age",
#'   padj = "fdr")
#' head(edge_result)
#' if(any(edge_result[,5]<0.05)){
#'   show(edge_result[which(edge_result[,5]<0.05),1:5])
#' }
#' }
#'
#' @importFrom stats as.formula coefficients lm p.adjust p.adjust.methods
#' @importFrom parallel detectCores makeCluster clusterExport parSapply stopCluster
#' @export

edge_lm <- function(net,
                   nnodes,
                   idata,
                   mod,
                   diag = FALSE,
                   padj,
                   cores = NULL,
                   expList = NULL,
                   verbose = TRUE,
                   ...){

  # Check 'diag' and 'nnodes'
  if(!is.logical(diag)) stop("STOP: diagonal variable must be logical!")
  if(nnodes%%1 != 0) stop("STOP: number of nodes must be integer!!")

  # Edge positions
  tri_pos <- which(upper.tri(matrix(nrow = nnodes, ncol = nnodes), diag = diag),
                   arr.ind = T)

  # Check if everything is alright
  net_dim <- dim(net)
  # If 3D input, reshape to 2D
  if(length(net_dim)==3){
    # Check number of variables
    if(!all(net_dim[1:2]==nnodes)) stop("STOP: number of matrix variables does not match volume dimmensions")
    # Extract edges
    mx <- sapply(1:nrow(tri_pos), function(x) net[tri_pos[x,1],tri_pos[x,2],])
  } else if(length(net_dim)==2){
    # Check number of variables
    if(net_dim[2] != nrow(tri_pos)) stop("STOP: number of matrix variables does not match column numbers")
    # Extract edges
    mx <- net
  } else stop("STOP: network input has to be a 2D or 3D array!!")
  if(nrow(idata)!=nrow(mx)) stop("STOP: network volume dimension 3 and dataset must match!!")
  if(!is.element(padj,p.adjust.methods)) stop("STOP: padj has to be a valid p.adjust.method!!")
  if(!is.numeric(as.matrix(mx))) stop("STOP: network array is not numeric!!")


  # Apply example model edgewise
  lm_form <- as.formula(paste0("mx[,",1,"]",mod))
  fit_ex <- coefficients(summary(lm(lm_form, idata, ...)))
  # Create empty object to store F-statistic and p-values edgewise
  obsT_ncol <- (nrow(fit_ex)-1)*3
  obsT <- matrix(as.numeric(NA), nrow = nrow(tri_pos), ncol = obsT_ncol)

  # Compute observed stats
  if(verbose) cat("Computing observed stats")
  # Set (or not) parallelization
  if(is.null(cores)){
    # Apply example model edgewise
    obsT[,(1:obsT_ncol)[-seq(3,obsT_ncol,3)]] <- t(sapply(1:nrow(tri_pos),
    function(x){
      lm_form <- as.formula(paste0("mx[,",x,"]",mod))
      fit <- coefficients(summary(lm(lm_form, idata, ...)))
      c(t(fit[2:nrow(fit),3:4]))
    }))
  } else{
    if(cores%%1 != 0) stop("STOP: number of cores must be integer")
    #library(parallel)
    if(cores > detectCores()) stop("STOP: input number of cores is to high")
    # Set 'clusterExport' varlist
    vlist <- c("tri_pos", "mod", "idata", "mx")
    if(!is.null(expList)){
      if(!is.character(expList)) stop("STOP: input 'expList' is not a character string")
      vlist <- c(vlist,expList)
    }
    if(verbose) cat(paste0(" (setting ",cores," cluster cores)"))
    cl <- makeCluster(cores)

    # Set variables for local function from the main function environment
    clusterExport(cl=cl, varlist=vlist, envir=environment())
    # Apply example model edgewise
    obsT[,(1:obsT_ncol)[-seq(3,obsT_ncol,3)]] <- t(parSapply(cl, 1:nrow(tri_pos),
    function(x){
      lm_form <- as.formula(paste0("mx[,",x,"]",mod))
      fit <- coefficients(summary(lm(lm_form, idata, ...)))
      c(t(fit[2:nrow(fit),3:4]))
    }))
    stopCluster(cl)
  }
  if(verbose) cat(".\n")

  # Applied multiple comparison correction with p.adjust
  if(verbose) cat("Computing p.adjust.method")
  obsT[,seq(3,obsT_ncol,3)] <- sapply(seq(2,obsT_ncol,3),function(x) p.adjust(p = obsT[,x], method = padj))
  # Set output list
  colnames(obsT) <- paste0(rep(rownames(fit_ex)[2:nrow(fit_ex)], each = 3),
                           c("_t","_p",paste0("_p",padj)))
  if(verbose) cat(".\n")
  
  # Add edge labels
  obsT <- cbind(tri_pos,obsT)
  colnames(obsT)[1:2] <- c("2Drow","2Dcol")
  
  # Return results
  return(as.data.frame(obsT))
}
