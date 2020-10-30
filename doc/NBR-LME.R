## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(NBR)
data("voles")
brain_labs <- NBR:::voles_roi
dim(voles)
head(voles)[1:8]

## ----upper triangle-----------------------------------------------------------
nnodes <- length(brain_labs)
tri_pos <- which(upper.tri(matrix(nrow = nnodes, ncol = nnodes)), arr.ind = T)
head(tri_pos)

## ----input network, fig.align = "center"--------------------------------------
library(lattice)
avg_mx <- matrix(0, nrow = nnodes, ncol = nnodes)
avg_mx[upper.tri(avg_mx)] <- apply(voles[-(1:3)], 2, function(x) mean(x, na.rm=TRUE))
avg_mx <- avg_mx + t(avg_mx)
# Set max-absolute value in order to set a color range centered in zero.
flim <- max(abs(avg_mx))
levelplot(avg_mx, main = "Average", ylab = "ROI", xlab = "ROI",
          at = seq(-flim, flim, length.out = 100))

## ----NBR-LME, eval = FALSE----------------------------------------------------
#  set.seed(18900217)
#  before <- Sys.time()
#  library(nlme)
#  nbr_result <- nbr_lme_aov(net = voles[,-(1:3)],
#    nnodes = 16,
#    idata = voles[,1:3],
#    nperm = 5,
#    mod = "~ Session*Sex",
#    rdm = "~ 1+Session|id",
#    na.action = na.exclude)
#  after <- Sys.time()
#  show(after-before)

## ----multicore NBR-LME, eval = FALSE------------------------------------------
#  set.seed(18900217)
#  before <- Sys.time()
#  library(nlme)
#  library(parallel)
#  nbr_result <- nbr_lme_aov(
#    net = voles[,-(1:3)],
#    nnodes = 16,
#    idata = voles[,1:3],
#    nperm = 1000,
#    nudist = T,
#    mod = "~ Session*Sex",
#    rdm = "~ 1+Session|id",
#    cores = detectCores(),
#    na.action = na.exclude
#    )
#  after <- Sys.time()
#  show(after-before)

## ----NBR-LME results----------------------------------------------------------
nbr_result <- NBR:::voles_nbr
show(nbr_result$fwe)

## ----component display, fig.align = "center"----------------------------------
# Plot significant edges
edge_mat <- array(0, dim(avg_mx))
edge_mat[nbr_result$components$Session[,2:3]] <- 1
levelplot(edge_mat, col.regions = rev(heat.colors(100)),
          main = "Component", ylab = "ROI", xlab = "ROI")

## ----component cum-pval, fig.height = 3, fig.width = 5, fig.align = "center"----
null_ses_str <- nbr_result$nudist[,2]  # Null distribution for Session strength
obs_ses_str <- nbr_result$fwe$Session[,4] # Observed Session strength
nperm <- length(null_ses_str)
cumpval <- cumsum(null_ses_str >= obs_ses_str)/(1:nperm)
# Plot p-value stability
plot(cumpval, type="l", ylim = c(0,0.06), las = 1,
           xlab = "Permutation index", ylab = "p-value",
           main = "Cumulative p-value for Session strength")
      abline(h=0.05, col="red", lty=2)
# Add binomial marginal error
mepval <- 2*sqrt(cumpval*(1-cumpval)/1:nperm)
lines(cumpval+mepval, col = "chartreuse4")
lines(cumpval-mepval, col = "chartreuse4")

