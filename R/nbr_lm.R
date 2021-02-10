#' @title Network-based R-statistics using Linear Model
#'
#' @description This function computes the specified linear model (LM) for each edge in the network,
#'  and calculates the family wise error (FWE) p-value for the size of the clusters of connected
#'  edges that are individually below the P threshold (\emph{thrP}), or above the T threshold
#'  (\emph{thrT}). FWE estimation is based on the null distribution of the maximum size of sets
#'  of connected edges (defined as above), obtained with \emph{nperm} permutations of the
#'  original data.
#'
#' @usage nbr_lm(net, nnodes, idata, mod, 
#'        alternative = c("two.sided", "less", "greater"),
#'        diag = FALSE, nperm, thrP = 0.05, thrT = NULL,
#'        cores = NULL, nudist = FALSE, expList = NULL, verbose = TRUE,
#'        ...)
#'
#' @param net 3D volume (2D matrices for each observation) or 2D matrix of edges as columns.
#' @param nnodes Number of network nodes.
#' @param idata Matrix or data.frame including independent variables of interest of the model.
#' @param mod Model, specify as a string, e.g., "~Group + Age".
#' @param alternative a character string specifying the alternative hypothesis, must be one of "two.sided" (default), "greater" or "less". You can specify just the initial letter.
#' @param diag Logical indicating if matrix diagonal is to be included in the analysis (default: FALSE).
#' @param nperm Number of permutations.
#' @param thrP Individual edge p-value threshold (if NULL, thrT should be given).
#' @param thrT Individual edge T-value threshold (if NULL, thrP should be given).
#' @param cores Number of selected cores for parallel computing (default: NULL).
#' @param nudist Logical indicating if null distribution should be returned (default: FALSE).
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
#' Regarding \emph{nperm}, I suggest first setting it to small values (5 or 10) in order to test that
#' everything runs fine. After that, set \emph{nperm} to 1000 or larger number to decrease the
#' margin of error of the FWE p-value (see \url{https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/Randomise/Theory#Conditional_Monte_Carlo_Permutation_Tests}
#' to explore the behavior of FWE p-value as a function of \emph{nperm}).
#'
#' @return List containing the observed statistics and their corresponding FWE p-values,
#'  if requested by \emph{nudist} it will return the null distribution.
#'  \enumerate{
#'   \item Observed statistics for every individual edge: corresponding subset of connected
#'    nodes and strength for each model term.
#'   \item FWE for components: binary and strength sum, with their corresponding FWE p-value.
#'   \item Null Distribution: maximal component size and strength for each permutation. Only
#'    returned if \emph{nudist} is TRUE.
#' }
#'
#' @examples
#' data(frontal2D)
#' \donttest{
#' nbr_result <- nbr_lm(net = frontal2D[,-(1:3)], nnodes = 28,
#'   idata = frontal2D[,1:3], mod = "~ Group + Sex * Age",
#'   thrP = NULL, thrT = 4, nperm = 5)
#' show(nbr_result)
#' }
#'
#' @importFrom stats as.formula coefficients lm qt aggregate
#' @importFrom parallel detectCores makeCluster clusterExport parSapply stopCluster
#' @export

nbr_lm <- function(net,
                   nnodes,
                   idata,
                   mod,
                   alternative = c("two.sided","less", "greater"),
                   diag = FALSE,
                   nperm,
                   thrP = 0.05,
                   thrT = NULL,
                   cores = NULL,
                   nudist = FALSE,
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
  if(nperm%%1 != 0) stop("STOP: number of permutations must be integer!!")
  if(!is.numeric(as.matrix(mx))) stop("STOP: network array is not numeric!!")
  alt <- match.arg(alternative, c("two.sided","less", "greater"))

  # Apply example model edgewise
  lm_form <- as.formula(paste0("mx[,",1,"]",mod))
  fit_ex <- coefficients(summary(lm(lm_form, idata, ...)))
  # Create empty object to store F-statistic and p-values edgewise
  obsT <- matrix(0, nrow = nrow(tri_pos), ncol = (nrow(fit_ex)-1)*2)

  # Compute observed stats
  if(verbose) cat("Computing observed stats")
  # Set (or not) parallelization
  if(is.null(cores)){
    # Apply example model edgewise
    obsT <- t(sapply(1:nrow(tri_pos), function(x){
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
    obsT <- t(parSapply(cl, 1:nrow(tri_pos), function(x){
      lm_form <- as.formula(paste0("mx[,",x,"]",mod))
      fit <- coefficients(summary(lm(lm_form, idata, ...)))
      c(t(fit[2:nrow(fit),3:4]))
    }))
    stopCluster(cl)
  }
  colnames(obsT) <- paste0(rep(rownames(fit_ex)[2:nrow(fit_ex)], each = 2),c("_t","_p"))
  if(verbose) cat(".\n")

  # Find components based on T-statistic or p-values thresholds
  if(all(is.null(thrT),is.null(thrP))) stop("STOP: t- or p-value threshold needed!")
  if(is.null(thrP)){thr <- thrT; thr_idx <- 1}
  if(is.null(thrT)){thr <- thrP; thr_idx <- 2}

  # Initial variables for component search
  obs_comp <- matrix(1:nnodes, nrow = nnodes, ncol = nrow(fit_ex)-1)
  obs_list <- vector("list",nrow(fit_ex)-1)
  # Set names
  names(obs_list) <- colnames(obs_comp) <- rownames(fit_ex)[2:nrow(fit_ex)]
  # Label component for each node and edge (and store strength as well)
  if(thr_idx == 1){
    # For each independent variable in the LM
    for(ii in 1:length(obs_list)){
      # Find t-statistics (t > thrT) edges
      if(alt=="two.sided"){
        edges <- which(abs(obsT[,ii*2-thr_idx])>thr)
      } else{
        if(alt=="less"){
          edges <- which(obsT[,ii*2-thr_idx]<thr)
        } else edges <- which(obsT[,ii*2-thr_idx]>thr)
      }
      if(length(edges)>0){
        # Generate empty object for component label
        component <- vector("integer", length(edges))
        # Store strength intensity
        strength <- abs(obsT[edges,ii*2-thr_idx])-abs(thr)
        # Mantain polarity
        strength <- strength*(1-(obsT[edges,ii*2-thr_idx]<0)*2)
        # Label nodes and edge component
        for(jj in 1:length(edges)) component[jj] <- obs_comp[tri_pos[edges[jj],],ii] <- min(obs_comp[tri_pos[edges[jj],],ii])
        # Store results
        obs_list[[ii]] <- cbind(edges,tri_pos[edges,1],tri_pos[edges,2],component,strength)
        colnames(obs_list[[ii]]) <- c("2Dcol","3Drow","3Dcol","comp","strn")
      }
    }
  } else{
    # For each independent variable in the LM
    for(ii in 1:length(obs_list)){
      # Find significant (p < thrP) edges
      if(alt=="two.sided"){
        edges <- which(obsT[,ii*thr_idx]<thr)
      } else{
        if(alt=="less"){
          edges <- which(obsT[,ii*thr_idx]<thr*2 & obsT[,ii*thr_idx-1]<0)
        } else edges <- which(obsT[,ii*thr_idx]<thr*2 & obsT[,ii*thr_idx-1]>0)
      }
      if(length(edges)>0){
        # Generate empty object for component label
        component <- vector("integer", length(edges))
        # Store strength intensity (bi-sided)
        if(alt=="two.sided"){
          strength <- abs(obsT[edges,ii*thr_idx-1])-qt(1-thrP/2,nrow(idata)-nrow(fit_ex))
        } else{
          strength <- abs(obsT[edges,ii*thr_idx-1])-qt(1-thrP,nrow(idata)-nrow(fit_ex))
        }
        # Mantain polarity
        strength <- strength*(1-(obsT[edges,ii*thr_idx-1]<0)*2)
        # Label nodes and edge component
        for(jj in 1:length(edges)){
          # Maximum label to be changed
          max_lab <- max(obs_comp[tri_pos[edges[jj],],ii])
          min_lab <- min(obs_comp[tri_pos[edges[jj],],ii])
          if(sum(component==max_lab)>0) component[which(component==max_lab)] <- min_lab
          pre_idx <- which(obs_comp[,ii]==max_lab)
          obs_comp[pre_idx,ii] <- component[jj] <- min_lab
        }
        # Store results
        obs_list[[ii]] <- cbind(edges,tri_pos[edges,1],tri_pos[edges,2],component,strength)
        colnames(obs_list[[ii]]) <- c("2Dcol","3Drow","3Dcol","comp","strn")
      }
    }
  }

  # Permutate dataset
  if(verbose) cat("Computing permutated stats")
  # Create empty object to store FWE null distribution
  null_dist <- array(0, dim=c(nperm, 2*length(obs_list)))
  colnames(null_dist) <- paste0(rep(names(obs_list), each=2),c("_ncomp","_strn"))
  # Set (or not) parallelization
  # Set (or not) parallelization
  if(!is.null(cores)){
    # Set 'clusterExport' varlist
    vlist <- c("tri_pos", "mod", "pdata", "mx")
    if(!is.null(expList)) vlist <- c(vlist,expList)
    if(verbose) cat(paste0(" (setting ",cores," cluster cores)"))
    cl <- makeCluster(cores)
  }
  if(verbose) cat(".\n")

  # Apply permutations
  if(verbose) cat("Permutation progress: ")
  for(pp in 1:nperm){

    # Print progress
    if(verbose) if(pp%%100 == 1) cat("...")
    if(verbose) if(pp%%100 == 0) cat(pp)

    # Create empty object to store F-statistic and p-values edgewise
    permT <- matrix(0, nrow = nrow(tri_pos), ncol = (nrow(fit_ex)-1)*2)
    # Permutate inference dataset
    pdata <- idata[sample(1:nrow(idata)),]

    # Set (or not) parallelization
    if(is.null(cores)){
      # Apply example model edgewise
      permT <- t(sapply(1:nrow(tri_pos), function(x){
        lm_form <- as.formula(paste0("mx[,",x,"]",mod))
        fit <- coefficients(summary(lm(lm_form, pdata, ...)))
        c(t(fit[2:nrow(fit),3:4]))
      }))
    } else{
      # Set variables for local function from the main function environment
      clusterExport(cl=cl, varlist=vlist, envir=environment())
      # Apply example model edgewise
      permT <- t(parSapply(cl, 1:nrow(tri_pos), function(x){
        lm_form <- as.formula(paste0("mx[,",x,"]",mod))
        fit <- coefficients(summary(lm(lm_form, pdata, ...)))
        c(t(fit[2:nrow(fit),3:4]))
      }))
    }

    # Initial variables for component search
    perm_comp <- matrix(1:nnodes, nrow = nnodes, ncol = nrow(fit_ex)-1)
    perm_list <- vector("list",nrow(fit_ex)-1)
    # Label component for each node and edge (and store strength as well)
    if(thr_idx == 1){
      # For each independent variable in the LM
      for(ii in 1:length(perm_list)){
        # Find t-statistics (t > thrT) edges
        if(alt=="two.sided"){
          edges <- which(abs(permT[,ii*2-thr_idx])>thr)
        } else{
          if(alt=="less"){
            edges <- which(permT[,ii*2-thr_idx]<thr)
          } else edges <- which(permT[,ii*2-thr_idx]>thr)
        }
        if(length(edges)>0){
          # Generate empty object for component label
          component <- vector("integer", length(edges))
          # Store strength
          strength <- abs(permT[edges,ii*2-thr_idx])-abs(thrT)
          # Label nodes and edge component
          for(jj in 1:length(edges)) component[jj] <- perm_comp[tri_pos[edges[jj],],ii] <- min(perm_comp[tri_pos[edges[jj],],ii])
          # Store results
          perm_list[[ii]] <- cbind(component,strength)
        }
      }
    } else{
      # For each independent variable in the LM
      for(ii in 1:length(perm_list)){
        # Find significant (p < thrP) edges
        if(alt=="two.sided"){
          edges <- which(permT[,ii*thr_idx]<thr)
        } else{
          if(alt=="less"){
            edges <- which(permT[,ii*thr_idx]<thr*2 & permT[,ii*thr_idx-1]<0)
          } else edges <- which(permT[,ii*thr_idx]<thr*2 & permT[,ii*thr_idx-1]>0)
        }
        if(length(edges)>0){
          # Generate empty object for component label
          component <- vector("integer", length(edges))
          # Store strength
          if(alt=="two.sided"){
            strength <- abs(permT[edges,ii*thr_idx-1])-qt(1-thrP/2,nrow(idata)-nrow(fit_ex))
          } else{
            strength <- abs(permT[edges,ii*thr_idx-1])-qt(1-thrP,nrow(idata)-nrow(fit_ex))
          }
          # Label nodes and edge component
          for(jj in 1:length(edges)) component[jj] <- perm_comp[tri_pos[edges[jj],],ii] <- min(perm_comp[tri_pos[edges[jj],],ii])
          # Store results
          perm_list[[ii]] <- cbind(component,strength)
        }
      }
    }

    # Save maximum values for the FWE null distribution
    for(ii in 1:length(perm_list)){
      if(!is.null(perm_list[[ii]])){
        null_dist[pp,(2*ii)-1] <- max(table(perm_list[[ii]][,1]))
        null_dist[pp,(2*ii)] <- max(aggregate(perm_list[[ii]][,2]~perm_list[[ii]][,1], perm_list[[ii]], sum)[,2])
      }
    }

  }#for(pp in 1:nperm)
  if(!is.null(cores)) stopCluster(cl)
  if(verbose) cat(".\n")

  # Compute FWE for each observed component
  # Quantile observed values within the null distribution
  fwe_list <- vector("list",nrow(fit_ex)-1)
  # Set names
  names(fwe_list) <- rownames(fit_ex)[2:nrow(fit_ex)]
  # Compute quantiles
  for(ii in 1:length(obs_list)){
    if(!is.null(obs_list[[ii]])){
      # Number of components FWE
      obs_ncomp <- ncompFWE <- table(obs_list[[ii]][,4])
      for(cc in 1:length(obs_ncomp)) ncompFWE[cc] <- sum(null_dist[,(2*ii)-1] > obs_ncomp[cc])/nperm
      # Components strength
      obs_strn <- strnFWE <- aggregate(obs_list[[ii]][,5]~obs_list[[ii]][,4], obs_list[[ii]], function(x) sum(abs(x)))
      for(cc in 1:nrow(obs_strn)) strnFWE[cc,2] <- sum(null_dist[,(2*ii)] > obs_strn[cc,2])/nperm
      # Save results
      DF <- as.data.frame(cbind(obs_strn[,1], c(obs_ncomp), c(ncompFWE), obs_strn[,2], strnFWE[,2]))
      names(DF) <- c("Component","ncomp","ncompFWE","strn","strnFWE")
      rownames(obs_strn) <- c()
      fwe_list[[ii]] <- DF
    }
  }

  # Return results
  if(nudist){
    return(list(components = obs_list, fwe = fwe_list, nudist = null_dist))
  } else{
    return(list(components = obs_list, fwe = fwe_list))
  }
}
