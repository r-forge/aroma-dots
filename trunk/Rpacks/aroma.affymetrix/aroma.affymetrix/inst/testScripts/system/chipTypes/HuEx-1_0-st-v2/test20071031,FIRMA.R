verbose <- Arguments$getVerbose(threshold=-3);
timestampOn(verbose);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setup data set
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
name <- "Affy_tissue_comm";
chipType <- "HuEx-1_0-st-v2";
coreChipType <- paste(chipType, "core", sep=",");
cdfCore <- AffymetrixCdfFile$fromChipType(coreChipType);

csTissue <- AffymetrixCelSet$fromName(name=name, chipType=chipType);
setCdf(csTissue, cdfCore);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Background correction and normalization
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
bcTissue <- RmaBackgroundCorrection(csTissue);
csBCTissue <- process(bcTissue,verbose=verbose);
qnTissue <- QuantileNormalization(csBCTissue, typesToUpdate="pm");
csNTissue <- process(qnTissue, verbose=verbose);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Probe-level summarization
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# if mergeGroups == FALSE, then each exon separately
plmTissue <- ExonRmaPlm(csNTissue, mergeGroups=TRUE);
print(plmTissue);

# if mergeGroups == FALSE, then each exon separately
plmNoMergeTissue <- ExonRmaPlm(csNTissue, mergeGroups=FALSE);
print(plmNoMergeTissue);

fit(plmTissue, verbose=verbose);
fit(plmNoMergeTissue, verbose=verbose);


# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Fit FIRMA
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
rsTissue <- calculateResiduals(plmTissue, verbose=verbose);
wsTissue <- calculateWeights(plmTissue,verbose=verbose);
rsNoMergeTissue <- calculateResiduals(plmNoMergeTissue, verbose=verbose);
wsNoMergeTissue <- calculateWeights(plmNoMergeTissue,verbose=verbose);

pafTissue <- getProbeAffinityFile(plmTissue);
cesTissue <- getChipEffectSet(plmTissue);
cesNoMergeTissue <- getChipEffectSet(plmNoMergeTissue);
pafNoMergeTissue <- getProbeAffinityFile(plmNoMergeTissue);

firmaTissue <- FirmaModel(plmTissue);
fit(firmaTissue, verbose=verbose);
fsTissue <- getFirmaScores(firmaTissue);
