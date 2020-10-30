#' @title Network-based R-statistics using Mixed Effects Models
#'
#' @description This function computes the specified (non)linear mixed models (LME) for each
#'  edge in the network, and calculates the family wise error (FWE) p-value for the size of the
#'  clusters of connected edges that are individually below the P threshold (\emph{thrP}), or
#'  above the T threshold (\emph{thrT}). FWE estimation is based on the null distribution of the
#'  maximum size of sets of connected edges (defined as above), obtained with \emph{nperm}
#'  permutations of the original data.
#'
#' @usage nbr_lme(net, nnodes, idata, mod, rdm, diag = FALSE,
#'         nperm, thrP = 0.05, thrT = NULL, cores = NULL,
#'         nudist = FALSE, expList = NULL, verbose = TRUE,
#'         ...)
#'
#' @param net 3D volume (2D matrices for each observation) or 2D matrix of edges as columns.
#' @param nnodes Number of network nodes.
#' @param idata Matrix or data.frame including independent variables of interest of the model.
#' @param mod Fixed effects, specify as a string, e.g., "~Session + Sex".
#' @param rdm Random effects, specify as a string, e.g., "~1+Session|id".
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
#' @inherit nbr_lm return details
#'
#' @examples
#' data(voles)
#' \donttest{
#' nbr_result <- nbr_lme(net = voles[,-(1:3)], nnodes = 16,
#'   idata = voles[,1:3], mod = "~ Session*Sex",
#'   rdm = "~ 1+Session|id", nperm = 5,
#'   na.action = na.exclude
#'   )
#' show(nbr_result)
#' }
#'
#' @importFrom stats as.formula coefficients qt ave aggregate
#' @importFrom parallel detectCores makeCluster clusterExport parSapply stopCluster
#' @importFrom nlme lme
#' @export

nbr_lme <- function(net,
                    nnodes,
                    idata,
                    mod,
                    rdm,
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
  idata_dim <- dim(idata)
  # If 3D input, reshape to 2D
  if(length(net_dim)==3){
    # Check number of variables
    if(!all(net_dim[1:2]==nnodes)) stop("STOP: number of matrix variables does not match volume dimmensions")
    # Extract edges
    mx <- sapply(1:nrow(tri_pos), function(x) net[tri_pos[x,1],tri_pos[x,2],])
    net_dim <- dim(mx)
  } else if(length(net_dim)==2){
    # Check number of variables
    if(net_dim[2] != nrow(tri_pos)) stop("STOP: number of matrix variables does not match column numbers")
    # Extract edges
    mx <- net
  } else stop("STOP: network input has to be a 2D or 3D array!!")
  if(net_dim[1]!=idata_dim[1]) stop("STOP: network volume dimension 3 and dataset must match!!")
  if(nperm%%1 != 0) stop("STOP: number of permutations must be integer!!")
  if(!is.numeric(as.matrix(mx))) stop("STOP: network array is not numeric!!")

  # Set input data.frame
  mx <- as.data.frame(mx)
  idata <- as.data.frame(cbind(idata,mx))
  # Apply example model edgewise
  ii <- 1
  lm_form <- as.formula(paste0(names(idata)[idata_dim[2]+ii], mod))
  fit_ex <- tryCatch(coefficients(summary(lme(lm_form, idata, as.formula(rdm), ...))),
                     error=function(e) NA)
  while(length(fit_ex)==1 & all(is.na(as.matrix(fit_ex)))){
    ii <- ii + 1
    if(ii > net_dim[2]) stop("STOP: fail to run every edgewise observed models.")
    lm_form <- as.formula(paste0(names(idata)[idata_dim[2]+ii], mod))
    fit_ex <- tryCatch(coefficients(summary(lme(lm_form, idata, as.formula(rdm), ...))),
                       error=function(e) NA)
  }
  # Create empty object to store F-statistic and p-values edgewise
  nedges <- nrow(tri_pos)
  nvars <- nrow(fit_ex)-1
  obsT <- matrix(0, nrow = nedges, ncol = nvars*2)

  # Compute observed stats
  if(verbose) cat("Computing observed stats")
  # Set (or not) parallelization
  if(is.null(cores)){
    # Apply example model edgewise
    obsT <- t(sapply(1:nedges, function(x){
      lm_form <- as.formula(paste0(names(idata)[idata_dim[2]+x], mod))
      fit <- tryCatch(coefficients(summary(lme(lm_form, idata, as.formula(rdm), ...))),
                      error=function(e) NA)
      # Store NA if individual model does not fit
      if(all(is.na(as.matrix(fit)))){
        return(rep(NA,nvars*2))
      } else{
        return(c(rbind(fit[1:nvars+1,4],fit[1:nvars+1,5])))
      }
    }))
  } else{
    if(cores%%1 != 0) stop("STOP: number of cores must be integer")
    #library(parallel)
    if(cores > detectCores()) stop("STOP: input number of cores is to high")
    # Set 'clusterExport' varlist
    vlist <- c("nedges", "nvars", "mod", "idata", "idata_dim", "lme", "rdm")
    if(!is.null(expList)){
      if(!is.character(expList)) stop("STOP: input 'expList' is not a character string")
      vlist <- c(vlist,expList)
    }
    if(verbose) cat(paste0(" (setting ",cores," cluster cores)"))
    cl <- makeCluster(cores)

    # Set variables for local function from the main function environment
    clusterExport(cl=cl, varlist=vlist, envir=environment())
    # Apply example model edgewise
    obsT <- t(parSapply(cl, 1:nedges, function(x){
      lm_form <- as.formula(paste0(names(idata)[idata_dim[2]+x], mod))
      fit <- tryCatch(coefficients(summary(lme(lm_form, idata, as.formula(rdm), ...))),
                      error=function(e) NA)
      # Store NA if individual model does not fit
      if(all(is.na(as.matrix(fit)))){
        return(rep(NA,nvars*2))
      } else{
        return(c(rbind(fit[1:nvars+1,4],fit[1:nvars+1,5])))
      }
    }))
    stopCluster(cl)
  }
  colnames(obsT) <- paste0(rep(rownames(fit_ex)[1:nvars+1], each = 2),c("_t","_p"))
  if(verbose) cat(".\n")

  # Find components based on F-statistic or p-values thresholds
  if(all(is.null(thrT),is.null(thrP))) stop("STOP: t- or p-value threshold needed!")
  if(is.null(thrP)){thr <- thrT; thr_idx <- 1}
  if(is.null(thrT)){thr <- thrP; thr_idx <- 2}

  # Initial variables for component search
  obs_comp <- matrix(1:nnodes, nrow = nnodes, ncol = nvars)
  obs_list <- vector("list",nvars)
  # Set names
  colnames(obs_comp) <- rownames(fit_ex)[1:nvars+1]
  names(obs_list) <- rownames(fit_ex)[1:nvars+1]
  # Label component for each node and edge (and store strength as well)
  if(thr_idx == 1){
    # For each independent variable in the LM
    for(ii in 1:nvars){
      # Find F-statistics (F > thrF) edges
      edges <- which(obsT[,ii*2-thr_idx]>thr)
      if(length(edges)>0){
        # Generate empty object for component label
        component <- vector("integer", length(edges))
        # Store strength intensity
        strength <- abs(obsT[edges,ii*2-thr_idx])-thr
        # Mantain polarity
        strength <- strength*(1-(obsT[edges,ii*2-thr_idx]<0)*2)
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
  } else{
    # For each independent variable in the LM
    for(ii in 1:nvars){
      # Find significant (p < thrP) edges
      edges <- which(obsT[,ii*thr_idx]<thr)
      if(length(edges)>0){
        # Generate empty object for component label
        component <- vector("integer", length(edges))
        # Store strength intensity (bi-sided)
        strength <- abs(obsT[edges,ii*thr_idx-1])-qt(1-thrP/2,fit_ex[ii+1,3])
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
  null_dist <- array(0, dim=c(nperm, 2*nvars))
  colnames(null_dist) <- paste0(rep(names(obs_list), each=2),c("_ncomp","_strn"))
  # Set (or not) parallelization
  if(!is.null(cores)){
    # Set 'clusterExport' varlist
    vlist <- c("nedges", "nvars", "mod", "pdata", "idata_dim", "lme", "rdm")
    if(!is.null(expList)) vlist <- c(vlist,expList)
    if(verbose) cat(paste0(" (setting ",cores," cluster cores)"))
    cl <- makeCluster(cores)
  }
  if(verbose) cat(".\n")

  # Find grouping variable
  rdmsplit <- unlist(strsplit(rdm,"[|]"))
  if(length(rdmsplit) != 2) stop("STOP: random formula must have one '|' character.")
  if(length(grep("/",rdmsplit[2]))>0) stop("STOP: this function does not allow nested random structures... yet.")
  # Set grouping structure index
  id_idx <- grep(rdmsplit[2],names(idata))

  # Apply permutations
  if(verbose) cat("Permutation progress: ")
  for(pp in 1:nperm){

    # Print progress
    if(verbose) if(pp%%100 == 1) cat("...")
    if(verbose) if(pp%%100 == 0) cat(pp)

    # Create empty object to store t-statistic and p-values edgewise
    permT <- matrix(0, nrow = nedges, ncol = nvars*2)
    # Permutate inference dataset (within block)
    perm_idx <- as.vector(sapply(as.data.frame(1:idata_dim[1]),
                                 function(x) ave(x, idata[,id_idx], FUN = sample)))
    pdata <- idata
    # Permute only original idata columns
    pdata[,1:idata_dim[2]] <- idata[perm_idx,1:idata_dim[2]]

    # Set (or not) parallelization
    if(is.null(cores)){
      # Apply example model edgewise
      permT <- t(sapply(1:nedges, function(x){
        lm_form <- as.formula(paste0(names(pdata)[idata_dim[2]+x], mod))
        fit <- tryCatch(coefficients(summary(lme(lm_form, pdata, as.formula(rdm), ...))),
                        error=function(e) NA)
        # Store NA if individual model does not fit
        if(all(is.na(as.matrix(fit)))){
          return(rep(NA,nvars*2))
        } else{
          return(c(rbind(fit[1:nvars+1,4],fit[1:nvars+1,5])))
        }
      }))
    } else{
      # Set variables for local function from the main function environment
      clusterExport(cl=cl, varlist=vlist, envir=environment())
      # Apply example model edgewise
      permT <- t(parSapply(cl, 1:nedges, function(x){
        lm_form <- as.formula(paste0(names(pdata)[idata_dim[2]+x], mod))
        fit <- tryCatch(coefficients(summary(lme(lm_form, pdata, as.formula(rdm), ...))),
                        error=function(e) NA)
        # Store NA if individual model does not fit
        if(all(is.na(as.matrix(fit)))){
          return(rep(NA,nvars*2))
        } else{
          return(c(rbind(fit[1:nvars+1,4],fit[1:nvars+1,5])))
        }
      }))
    }

    # Initial variables for component search
    perm_comp <- matrix(1:nnodes, nrow = nnodes, ncol = nvars)
    perm_list <- vector("list", nvars)
    # Label component for each node and edge (and store strength as well)
    if(thr_idx == 1){
      # For each independent variable in the LM
      for(ii in 1:nvars){
        # Find F-statistics (F > thrF) edges
        edges <- which(permT[,ii*2-thr_idx]>thr)
        if(length(edges)>0){
          # Generate empty object for component label
          component <- vector("integer", length(edges))
          # Store strength
          strength <- abs(permT[edges,ii*2-thr_idx])-thrT
          # Label nodes and edge component
          for(jj in 1:length(edges)){
            # Maximum label to be changed
            max_lab <- max(perm_comp[tri_pos[edges[jj],],ii])
            min_lab <- min(perm_comp[tri_pos[edges[jj],],ii])
            if(sum(component==max_lab)>0) component[which(component==max_lab)] <- min_lab
            pre_idx <- which(perm_comp[,ii]==max_lab)
            perm_comp[pre_idx,ii] <- component[jj] <- min_lab
          }
          # Store results
          perm_list[[ii]] <- cbind(component,strength)
        }
      }
    } else{
      # For each independent variable in the LM
      for(ii in 1:nvars){
        # Find significant (p < thrP) edges
        edges <- which(permT[,ii*thr_idx]<thr)
        if(length(edges)>0){
          # Generate empty object for component label
          component <- vector("integer", length(edges))
          # Store strength
          strength <- abs(permT[edges,ii*thr_idx-1])-qt(1-thrP/2,fit_ex[ii+1,3])
          # Label nodes and edge component
          for(jj in 1:length(edges)){
            # Maximum label to be changed
            max_lab <- max(perm_comp[tri_pos[edges[jj],],ii])
            min_lab <- min(perm_comp[tri_pos[edges[jj],],ii])
            if(sum(component==max_lab)>0) component[which(component==max_lab)] <- min_lab
            pre_idx <- which(perm_comp[,ii]==max_lab)
            perm_comp[pre_idx,ii] <- component[jj] <- min_lab
          }
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
  fwe_list <- vector("list", nvars)
  # Set names
  names(fwe_list) <- rownames(fit_ex)[1:nvars+1]
  # Compute quantiles
  for(ii in 1:nvars){
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
