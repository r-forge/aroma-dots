library(aroma.affymetrix);

log <- Arguments$getVerbose(-3);
timestampOn(log);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup data set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
name <- "Affymetrix-HeartBrain";
chipType <- "HuEx-1_0-st-v2";
coreChipType <- paste(chipType, "core", sep=",");
cdfCore <- AffymetrixCdfFile$fromChipType(coreChipType);

csTissue <- AffymetrixCelSet$fromName(name=name, chipType=chipType);
setCdf(csTissue, cdfCore);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Background correction and normalization
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
bcTissue <- RmaBackgroundCorrection(csTissue);
csBCTissue <- process(bcTissue,verbose=log);
qnTissue <- QuantileNormalization(csBCTissue, typesToUpdate="pm");
csNTissue <- process(qnTissue, verbose=log);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Probe-level summarization
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# if mergeGroups == FALSE, then each exon separately
plmTissue <- ExonRmaPlm(csNTissue, mergeGroups=TRUE);
print(plmTissue);

# if mergeGroups == FALSE, then each exon separately
plmNoMergeTissue <- ExonRmaPlm(csNTissue, mergeGroups=FALSE);
print(plmNoMergeTissue);

fit(plmTissue, verbose=log);
fit(plmNoMergeTissue, verbose=log);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Fit FIRMA
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
rsTissue <- calculateResiduals(plmTissue, verbose=log);
wsTissue <- calculateWeights(plmTissue,verbose=log);
rsNoMergeTissue <- calculateResiduals(plmNoMergeTissue, verbose=log);
wsNoMergeTissue <- calculateWeights(plmNoMergeTissue,verbose=log);

pafTissue <- getProbeAffinityFile(plmTissue);
cesTissue <- getChipEffectSet(plmTissue);
cesNoMergeTissue <- getChipEffectSet(plmNoMergeTissue);
pafNoMergeTissue <- getProbeAffinityFile(plmNoMergeTissue);

firmaTissue <- FirmaModel(plmTissue);
fit(firmaTissue, verbose=log);
fsTissue <- getFirmaScores(firmaTissue);
