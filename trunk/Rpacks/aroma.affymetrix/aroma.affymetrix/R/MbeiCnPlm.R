###########################################################################/**
# @RdocClass MbeiCnPlm
#
# @title "The MbeiCnPlm class"
#
# \description{
#  @classhierarchy
#
#  This class represents the log-additive model used in RMA.
#  It can be used to fit the model on a @see "AffymetrixCelSet".
# }
# 
# @synopsis
#
# \arguments{
#   \item{...}{Arguments passed to @see "MbeiSnpPlm".}
# }
#
# \section{Fields and Methods}{
#  @allmethods "public"
# }
#
# \section{Model}{
#   Consider a specific unit group.  The log-additive model in RMA is:
#
#    \deqn{log_2(y_{ij})} = \beta_i + \alpha_j + \eps_{ij}
#
#   where \eqn{\beta_i} are the chip effects for arrays \eqn{i=1,...,I}, 
#   and \eqn{\alpha_j} are the probe affinities for probes \eqn{j=1,...,J}.
#   The \eqn{\eps_{ij}} are zero-mean noise with equal variance.
#
#   To minimize the risk for mistakes, all probe-affinity models return
#   parameter estimates on the intensity scale.  That is, this class
#   returns \eqn{\theta_i = 2^\beta_i} and \eqn{\phi_i = 2^\alpha_i},
#   cf. the multiplicative model of Li & Wong.
#
#   Use @seemethod "getProbeAffinities" to get the probe-affinity estimates.
# }
#
# @author
#
# \section{Model estimates}{
#   The estimated probe affinities are represented by the
#   @see "CnProbeAffinityFile" class.  
# }
#
# \references{
# }
#*/###########################################################################
setConstructorS3("MbeiCnPlm", function(..., combineAlleles=FALSE, name="modelMbeiCnPlm") {
  extend(MbeiSnpPlm(..., name=name), c("MbeiCnPlm", uses(CnPlm())),
    combineAlleles = combineAlleles
  )
})


############################################################################
# HISTORY:
# 2006-09-12
# o Recreated.
############################################################################
