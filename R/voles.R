#' Prairie voles functional connectivity
#'
#' A dataset containing the functional connectivity between
#' 16 brain areas of 32 prairie voles in three different sessions.
#'
#' Data is based on an experiment of social bonding in prairie voles.
#' Functional connectivity was measured as the Pearson correlation
#' between the average fMRI signal from the regions of interest (ROI)
#' within 16 anatomical areas of brain. Then, a total of
#' 120 pairwise connections are contained in the dataset.
#' NOTE: This is not the original data of the study!
#'
#' @format A data.frame with 96 rows and 123 variables:
#' \describe{
#'   \item{id}{Subject ID, factor.}
#'   \item{Sex}{Factor: female (F) or male (M).}
#'   \item{Session}{Factor: 1st, 2nd, or 3rd.}
#'   \item{ACC.AON}{Functional connectivity between ACC and AON regions, numeric}
#'   \item{ACC.BLA}{Functional connectivity between ACC and BLA regions, numeric}
#'   \item{AON.BLA}{Functional connectivity between AON and BLA regions, numeric}
#'   \item{ACC.BNST}{Functional connectivity between ACC and BNST regions, numeric}
#'   \item{AON.BNST}{Functional connectivity between AON and BNST regions, numeric}
#'   \item{BLA.BNST}{Functional connectivity between BLA and BNST regions, numeric}
#'   \item{ACC.LS}{Functional connectivity between ACC and LS regions, numeric}
#'   \item{AON.LS}{Functional connectivity between AON and LS regions, numeric}
#'   \item{BLA.LS}{Functional connectivity between BLA and LS regions, numeric}
#'   \item{BNST.LS}{Functional connectivity between BNST and LS regions, numeric}
#'   \item{ACC.MeA}{Functional connectivity between ACC and MeA regions, numeric}
#'   \item{AON.MeA}{Functional connectivity between AON and MeA regions, numeric}
#'   \item{BLA.MeA}{Functional connectivity between BLA and MeA regions, numeric}
#'   \item{BNST.MeA}{Functional connectivity between BNST and MeA regions, numeric}
#'   \item{LS.MeA}{Functional connectivity between LS and MeA regions, numeric}
#'   \item{ACC.MOB}{Functional connectivity between ACC and MOB regions, numeric}
#'   \item{AON.MOB}{Functional connectivity between AON and MOB regions, numeric}
#'   \item{BLA.MOB}{Functional connectivity between BLA and MOB regions, numeric}
#'   \item{BNST.MOB}{Functional connectivity between BNST and MOB regions, numeric}
#'   \item{LS.MOB}{Functional connectivity between LS and MOB regions, numeric}
#'   \item{MeA.MOB}{Functional connectivity between MeA and MOB regions, numeric}
#'   \item{ACC.mPFC}{Functional connectivity between ACC and mPFC regions, numeric}
#'   \item{AON.mPFC}{Functional connectivity between AON and mPFC regions, numeric}
#'   \item{BLA.mPFC}{Functional connectivity between BLA and mPFC regions, numeric}
#'   \item{BNST.mPFC}{Functional connectivity between BNST and mPFC regions, numeric}
#'   \item{LS.mPFC}{Functional connectivity between LS and mPFC regions, numeric}
#'   \item{MeA.mPFC}{Functional connectivity between MeA and mPFC regions, numeric}
#'   \item{MOB.mPFC}{Functional connectivity between MOB and mPFC regions, numeric}
#'   \item{ACC.NAcc}{Functional connectivity between ACC and NAcc regions, numeric}
#'   \item{AON.NAcc}{Functional connectivity between AON and NAcc regions, numeric}
#'   \item{BLA.NAcc}{Functional connectivity between BLA and NAcc regions, numeric}
#'   \item{BNST.NAcc}{Functional connectivity between BNST and NAcc regions, numeric}
#'   \item{LS.NAcc}{Functional connectivity between LS and NAcc regions, numeric}
#'   \item{MeA.NAcc}{Functional connectivity between MeA and NAcc regions, numeric}
#'   \item{MOB.NAcc}{Functional connectivity between MOB and NAcc regions, numeric}
#'   \item{mPFC.NAcc}{Functional connectivity between mPFC and NAcc regions, numeric}
#'   \item{ACC.PVN}{Functional connectivity between ACC and PVN regions, numeric}
#'   \item{AON.PVN}{Functional connectivity between AON and PVN regions, numeric}
#'   \item{BLA.PVN}{Functional connectivity between BLA and PVN regions, numeric}
#'   \item{BNST.PVN}{Functional connectivity between BNST and PVN regions, numeric}
#'   \item{LS.PVN}{Functional connectivity between LS and PVN regions, numeric}
#'   \item{MeA.PVN}{Functional connectivity between MeA and PVN regions, numeric}
#'   \item{MOB.PVN}{Functional connectivity between MOB and PVN regions, numeric}
#'   \item{mPFC.PVN}{Functional connectivity between mPFC and PVN regions, numeric}
#'   \item{NAcc.PVN}{Functional connectivity between NAcc and PVN regions, numeric}
#'   \item{ACC.RSC}{Functional connectivity between ACC and RSC regions, numeric}
#'   \item{AON.RSC}{Functional connectivity between AON and RSC regions, numeric}
#'   \item{BLA.RSC}{Functional connectivity between BLA and RSC regions, numeric}
#'   \item{BNST.RSC}{Functional connectivity between BNST and RSC regions, numeric}
#'   \item{LS.RSC}{Functional connectivity between LS and RSC regions, numeric}
#'   \item{MeA.RSC}{Functional connectivity between MeA and RSC regions, numeric}
#'   \item{MOB.RSC}{Functional connectivity between MOB and RSC regions, numeric}
#'   \item{mPFC.RSC}{Functional connectivity between mPFC and RSC regions, numeric}
#'   \item{NAcc.RSC}{Functional connectivity between NAcc and RSC regions, numeric}
#'   \item{PVN.RSC}{Functional connectivity between PVN and RSC regions, numeric}
#'   \item{ACC.VP}{Functional connectivity between ACC and VP regions, numeric}
#'   \item{AON.VP}{Functional connectivity between AON and VP regions, numeric}
#'   \item{BLA.VP}{Functional connectivity between BLA and VP regions, numeric}
#'   \item{BNST.VP}{Functional connectivity between BNST and VP regions, numeric}
#'   \item{LS.VP}{Functional connectivity between LS and VP regions, numeric}
#'   \item{MeA.VP}{Functional connectivity between MeA and VP regions, numeric}
#'   \item{MOB.VP}{Functional connectivity between MOB and VP regions, numeric}
#'   \item{mPFC.VP}{Functional connectivity between mPFC and VP regions, numeric}
#'   \item{NAcc.VP}{Functional connectivity between NAcc and VP regions, numeric}
#'   \item{PVN.VP}{Functional connectivity between PVN and VP regions, numeric}
#'   \item{RSC.VP}{Functional connectivity between RSC and VP regions, numeric}
#'   \item{ACC.VTA}{Functional connectivity between ACC and VTA regions, numeric}
#'   \item{AON.VTA}{Functional connectivity between AON and VTA regions, numeric}
#'   \item{BLA.VTA}{Functional connectivity between BLA and VTA regions, numeric}
#'   \item{BNST.VTA}{Functional connectivity between BNST and VTA regions, numeric}
#'   \item{LS.VTA}{Functional connectivity between LS and VTA regions, numeric}
#'   \item{MeA.VTA}{Functional connectivity between MeA and VTA regions, numeric}
#'   \item{MOB.VTA}{Functional connectivity between MOB and VTA regions, numeric}
#'   \item{mPFC.VTA}{Functional connectivity between mPFC and VTA regions, numeric}
#'   \item{NAcc.VTA}{Functional connectivity between NAcc and VTA regions, numeric}
#'   \item{PVN.VTA}{Functional connectivity between PVN and VTA regions, numeric}
#'   \item{RSC.VTA}{Functional connectivity between RSC and VTA regions, numeric}
#'   \item{VP.VTA}{Functional connectivity between VP and VTA regions, numeric}
#'   \item{ACC.Dent}{Functional connectivity between ACC and Dent regions, numeric}
#'   \item{AON.Dent}{Functional connectivity between AON and Dent regions, numeric}
#'   \item{BLA.Dent}{Functional connectivity between BLA and Dent regions, numeric}
#'   \item{BNST.Dent}{Functional connectivity between BNST and Dent regions, numeric}
#'   \item{LS.Dent}{Functional connectivity between LS and Dent regions, numeric}
#'   \item{MeA.Dent}{Functional connectivity between MeA and Dent regions, numeric}
#'   \item{MOB.Dent}{Functional connectivity between MOB and Dent regions, numeric}
#'   \item{mPFC.Dent}{Functional connectivity between mPFC and Dent regions, numeric}
#'   \item{NAcc.Dent}{Functional connectivity between NAcc and Dent regions, numeric}
#'   \item{PVN.Dent}{Functional connectivity between PVN and Dent regions, numeric}
#'   \item{RSC.Dent}{Functional connectivity between RSC and Dent regions, numeric}
#'   \item{VP.Dent}{Functional connectivity between VP and Dent regions, numeric}
#'   \item{VTA.Dent}{Functional connectivity between VTA and Dent regions, numeric}
#'   \item{ACC.HipD}{Functional connectivity between ACC and HipD regions, numeric}
#'   \item{AON.HipD}{Functional connectivity between AON and HipD regions, numeric}
#'   \item{BLA.HipD}{Functional connectivity between BLA and HipD regions, numeric}
#'   \item{BNST.HipD}{Functional connectivity between BNST and HipD regions, numeric}
#'   \item{LS.HipD}{Functional connectivity between LS and HipD regions, numeric}
#'   \item{MeA.HipD}{Functional connectivity between MeA and HipD regions, numeric}
#'   \item{MOB.HipD}{Functional connectivity between MOB and HipD regions, numeric}
#'   \item{mPFC.HipD}{Functional connectivity between mPFC and HipD regions, numeric}
#'   \item{NAcc.HipD}{Functional connectivity between NAcc and HipD regions, numeric}
#'   \item{PVN.HipD}{Functional connectivity between PVN and HipD regions, numeric}
#'   \item{RSC.HipD}{Functional connectivity between RSC and HipD regions, numeric}
#'   \item{VP.HipD}{Functional connectivity between VP and HipD regions, numeric}
#'   \item{VTA.HipD}{Functional connectivity between VTA and HipD regions, numeric}
#'   \item{Dent.HipD}{Functional connectivity between Dent and HipD regions, numeric}
#'   \item{ACC.HipV}{Functional connectivity between ACC and HipV regions, numeric}
#'   \item{AON.HipV}{Functional connectivity between AON and HipV regions, numeric}
#'   \item{BLA.HipV}{Functional connectivity between BLA and HipV regions, numeric}
#'   \item{BNST.HipV}{Functional connectivity between BNST and HipV regions, numeric}
#'   \item{LS.HipV}{Functional connectivity between LS and HipV regions, numeric}
#'   \item{MeA.HipV}{Functional connectivity between MeA and HipV regions, numeric}
#'   \item{MOB.HipV}{Functional connectivity between MOB and HipV regions, numeric}
#'   \item{mPFC.HipV}{Functional connectivity between mPFC and HipV regions, numeric}
#'   \item{NAcc.HipV}{Functional connectivity between NAcc and HipV regions, numeric}
#'   \item{PVN.HipV}{Functional connectivity between PVN and HipV regions, numeric}
#'   \item{RSC.HipV}{Functional connectivity between RSC and HipV regions, numeric}
#'   \item{VP.HipV}{Functional connectivity between VP and HipV regions, numeric}
#'   \item{VTA.HipV}{Functional connectivity between VTA and HipV regions, numeric}
#'   \item{Dent.HipV}{Functional connectivity between Dent and HipV regions, numeric}
#'   \item{HipD.HipV}{Functional connectivity between HipD and HipV regions, numeric}
#' }
#' @source \url{https://www.biorxiv.org/content/10.1101/752345v2}
"voles"
