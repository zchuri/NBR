## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup--------------------------------------------------------------------
library(NBR)
cmx <- NBR:::frontal3D          # Load 3D array
brain_labs <- NBR:::frontal_roi # Load node labels
phen <- NBR:::frontal_phen      # Load phenotypic info
dim(cmx)                        # Show 3D array dimensions

## ----input networks, fig.align = "center"-------------------------------------
library(lattice)
avg_mx <- apply(cmx, 1:2, mean)
# Set max-absolute value in order to set a color range centered in zero.
flim <- max(abs(avg_mx)[is.finite(avg_mx)])
levelplot(avg_mx, main = "Average", ylab = "ROI", xlab = "ROI",
          at = seq(-flim, flim, length.out = 100))

## ----input phenotypic info----------------------------------------------------
head(phen)
nrow(phen)
identical(nrow(phen), dim(cmx)[3])

## ----group-based NBR----------------------------------------------------------
set.seed(18900217) # Because R. Fisher is my hero
before <- Sys.time()
nbr_group <- nbr_lm_aov(net = cmx, nnodes = 28, idata = phen,
   mod = "~ Group", thrP = 0.01, nperm = 10)
after <- Sys.time()
show(after-before)

## ----multicore group-based NBR, eval = FALSE----------------------------------
#  set.seed(18900217)
#  library(parallel)
#  before <- Sys.time()
#  nbr_group <- nbr_lm_aov(net = cmx, nnodes = 28, idata = phen,
#     mod = "~ Group", thrP = 0.01, nperm = 100, cores = detectCores())
#  after <- Sys.time()
#  length(nbr_group)

## ----component display, fig.align = "center"----------------------------------
# Plot significant component
edge_mat <- array(0, dim(avg_mx))
edge_mat[nbr_group$components$Group[,2:3]] <- 1
levelplot(edge_mat, col.regions = rev(heat.colors(100)),
          main = "Component", ylab = "ROI", xlab = "ROI")
show(nbr_group$fwe$Group)

