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
print(cdfCore);

csTissue <- AffymetrixCelSet$fromName(name=name, chipType=chipType);
setCdf(csTissue, cdfCore);
print(csTissue);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Background correction and normalization
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
bc <- RmaBackgroundCorrection(csTissue);
print(bc);
csBCTissue <- process(bc, verbose=log);
print(csBCTissue);

qn <- QuantileNormalization(csBCTissue, typesToUpdate="pm");
print(qn);
csNTissue <- process(qn, verbose=log);
print(csNTissue);


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
# FIRMA
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
rsTissue <- calculateResiduals(plmTissue, verbose=log);
wsTissue <- calculateWeights(plmTissue,verbose=log);
rsNoMergeTissue <- calculateResiduals(plmNoMergeTissue, verbose=log);
wsNoMergeTissue <- calculateWeights(plmNoMergeTissue,verbose=log);

pafTissue <- getProbeAffinityFile(plmTissue);
cesTissue <- getChipEffectSet(plmTissue);
cesNoMergeTissue <- getChipEffectSet(plmNoMergeTissue);
pafNoMergeTissue <- getProbeAffinityFile(plmNoMergeTissue);

firma <- FirmaModel(plmTissue);
print(firma);

fit(firma, verbose=log);

fsTissue <- getFirmaScores(firma);
print(fsTissue);
