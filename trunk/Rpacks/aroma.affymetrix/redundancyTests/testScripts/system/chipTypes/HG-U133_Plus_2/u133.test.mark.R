
# --------------------
# need these files and directory structure ... download from (tissues):
# http://www.affymetrix.com/support/technical/sample_data/exon_array_data.affx
# http://www.affymetrix.com/support/technical/sample_data/hugene_1_0_array_data.affx
# --------------------
#-rw-r-----  1 mrobinson lab0605 13555904 Jun 19 14:03 rawData/heart_brain/HG-U133_Plus_2/u1332plus_ivt_cerebellum_A.CEL
#-rw-r-----  1 mrobinson lab0605 13550687 Jun 19 14:03 rawData/heart_brain/HG-U133_Plus_2/u1332plus_ivt_cerebellum_B.CEL
#-rw-r-----  1 mrobinson lab0605 13551860 Jun 19 14:03 rawData/heart_brain/HG-U133_Plus_2/u1332plus_ivt_cerebellum_C.CEL
#-rw-r-----  1 mrobinson lab0605 13554731 Jun 19 14:03 rawData/heart_brain/HG-U133_Plus_2/u1332plus_ivt_heart_A.CEL
#-rw-r-----  1 mrobinson lab0605 13553255 Jun 19 14:03 rawData/heart_brain/HG-U133_Plus_2/u1332plus_ivt_heart_B.CEL
#-rw-r-----  1 mrobinson lab0605 13551203 Jun 19 14:03 rawData/heart_brain/HG-U133_Plus_2/u1332plus_ivt_heart_C.CEL
#-rw-r--r--  1 mrobinson lab0605 66429204 Oct 11  2006 rawData/heart_brain/HuEx-1_0-st-v2/huex_wta_cerebellum_A.CEL
#-rw-r--r--  1 mrobinson lab0605 66407841 Oct 11  2006 rawData/heart_brain/HuEx-1_0-st-v2/huex_wta_cerebellum_B.CEL
#-rw-r--r--  1 mrobinson lab0605 66194497 Oct 11  2006 rawData/heart_brain/HuEx-1_0-st-v2/huex_wta_cerebellum_C.CEL
#-rw-r--r--  1 mrobinson lab0605 66595752 Oct 11  2006 rawData/heart_brain/HuEx-1_0-st-v2/huex_wta_heart_A.CEL
#-rw-r--r--  1 mrobinson lab0605 66360744 Oct 11  2006 rawData/heart_brain/HuEx-1_0-st-v2/huex_wta_heart_B.CEL
#-rw-r--r--  1 mrobinson lab0605 66797339 Oct 11  2006 rawData/heart_brain/HuEx-1_0-st-v2/huex_wta_heart_C.CEL

coef.tol<-.001

library(aroma.affymetrix)
getExprs <- function(ces, ...) {
  cdf <- getCdf(ces);
  theta <- extractMatrix(ces, ..., returnUgcMap=TRUE);
  ugcMap <- attr(theta, "unitGroupCellMap");
  rownames(theta) <- getUnitNames(cdf, ugcMap[,"unit"]);
  log2(theta)
}

# --------------------
# run aroma.affymetrix 
# --------------------
verbose <- Arguments$getVerbose(-8); timestampOn(verbose)
name <- "heart_brain"
chip1<- "HG-U133_Plus_2"
cs1 <- AffymetrixCelSet$fromName(name, chipType=chip1)
bc1 <- RmaBackgroundCorrection(cs1)
csBC1 <- process(bc1, verbose=verbose)
qn1 <- QuantileNormalization(csBC1, typesToUpdate="pm")
csN1 <- process(qn1, verbose=verbose) #time required
plm1.merge <- RmaPlm(csN1)
fit(plm1.merge, verbose=verbose) #time required
chp1.merge<-getChipEffectSet(plm1.merge)
d1<-getExprs(chp1.merge)

# --------------------
# run affyPLM to match
# --------------------
fn.u133 <- getPathnames(cs1)
raw<-ReadAffy(filenames=fn.u133)
plm<-fitPLM(raw,verbos=9)
d<-coefs(plm)

#png("rma.vs.RmaPlm.png")
#plot(d[,1],d1[m,1],pch=19,xlab="fitPLM-affyPLM",ylab="RmaBackgroundCorrection+RmaPlm-aroma.affymetrix")
#abline(0,1,col="blue")
#dev.off()

# --------------------
# check chip effects
# --------------------
m<-match(rownames(d),rownames(d1))
e<-d-d1[m,]
stopifnot(mean(e^2)<coef.tol)


# --------------------
# do more tests
# --------------------
