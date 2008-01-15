#########################################################################/**
# @RdocDocumentation "Weights"
#
# \description{
#   For several normalization and calibration methods the estimation of the
#   normalization (or calibration) function can be done with weights. 
#   Commonly, weights are proportional to a quality measure, that is,
#   the less quality we assign to a signal, the less influence (weight) it 
#   should have on the estimation of calibration and normalization 
#   functions.
# }
#
#
# \section{General about weights}{
#   The definition of a \emph{weight} is a single value in [0,1].
#   Weights outside this range and @NAs (missing values) are not allowed.
#
#   Below, we will define different entities such as signals, probes/spots, 
#   probesets, channels, arrays, and data-points.
#   To any of these entities weights may be assigned.
# }
#
#
# \section{Signals and signal weights}{
#   A \emph{signal} is a single value. 
#   A \emph{signal weight} is a weight assigned to a signal.
#   Thus, it is for entities within an array and never across/between arrays.
#
#   \emph{Example}: In two-color microarray data, there are two signals
#   for each spot, i.e. the red or the green signals, and each of them
#   can be assigned a different signal weight.
#   Typically, such signal weights are represented by an Nx2 @matrix,
#   where N is the number of probes/spots on the array.
#
#   \emph{Example}: In Affymetrix microarray data, which is single-channel
#   data, there is one signal per probe (in turn part of a probe set). 
#   Each such probe can be assigned a signal weight.
#   Typically, such signal weights are represented by an Nx1 @matrix,
#   where N is the number of probes/spots on the array.
# }
#
#
# \section{Probes and probe weights}{
#   A \emph{probe} is the smallest entity (not considering image pixels) 
#   on the array that measures the amount of hybridized samples in
#   one or several channels.
#
#   \emph{Example}: For two-color microarrays, a probe is a spot.
#   \emph{Example}: For Affymetrix arrays, a probe can be either a
#   perfect match probe (PM) or a mismatch probe (MM).
#
#   A \emph{probe weight} is a weight assigned to a probe/spot 
#   (not a probe set). 
#
#   \emph{Example}: For two-color data, the signals in the two channels 
#   for a given spot share the same probe weight.
#
#   \emph{Example}: For four-color data, the signals in the four channels 
#   for a given spot share the same probe weight.
#
#   \emph{Example}: For single-channel data such as Affymetrix data, the
#   probe weight is identical to a signal weight.
#
#   Typically, above signal weights are represented by an Nx1 @matrix,
#   where N is the number of probes/spots on the array.
#
#   The probe weight of probe $i$ must be equal to the mean of its
#   signal weights.
# }
#
#
# \section{Probesets and probeset weights}{
#   A \emph{probeset} consists of a set of probes. 
#
#   \emph{Example}: For two-color microarrays, probesets are not defined.
#   \emph{Example}: For Affymetrix arrays, a probeset is the set of
#   perfect match (PM) and mismatch (MM) probes corresponding to the
#   same gene. 
#
#   A \emph{probeset weight} is a weight assigned to a probeset.
#
#   Since Affymetrix is single-channel arrays, typically the above probeset 
#   weights are represented by an Nx1 @matrix, where N is the number of 
#   probesets.
#
#   The probeset weight for probeset $j$ must be equal to the mean of
#   its probe weights. (==signal weights)
#   by averaging the probe weights for each probeset.
# }
#
#
# \section{Data points and data-point weights}{
#   The definition of a \emph{data point} depends on the context.
#   It may be assigned to entities within an array, but also across/between
#   arrays.
#
#   \emph{Example}: (Paired data-point weight). In paired-channel 
#   normalization, such as curve-fit normalization (a.k.a. lowess intensity
#   normalization), two and only two channels are normalized together at
#   the same time, e.g. red and the green channels in two-color data, or two
#   two single-channel data set obtained from two different Affymetrix
#   arrays. Here a data point is constituted by two signals, e.g. 
#   \eqn{X_i = (G_i,R_i)}.
#   A data-point weight is assigned to the pair of signals corresponding to
#   the same spot or gene, e.g. for lowess normalization a data-point
#   weight is assigned to a log-ratio and a log-intensity.
#
#   \emph{Example}: (Multi-channel data point weight). In multi-channel
#   normalization, such as affine normalization or quantile normalization
#   (within a singel array and/or across multiple arrays), each data point
#   is consituted by multiple signals, e.g. for K two-color arrays it is
#   \eqn{X_i = (R[i,1],G[i,1],...,R[i,K],G[i,K])}.
#   To this data point, a \emph{data-point weight} can be assigned.
#
#   Typically, above signal weights are represented by an Nx1 @matrix,
#   where N is the number of data points.
#
#   Data-points weights can be generated from signal or probe weights, by
#   averaging them for each data point.
#
#   If not stated elsewise, arguments named \code{weights} are assumed
#   to take data-point weights.
# }
#
#
# \section{Arrays and array weights}{
#   An \emph{array weight} is a weight assigned to an array, that is,
#   to the complete set of signals in all channels constituting an array.
#
#   By definition, a channel weight can never apply across/between arrays.
#
#   \emph{Constraints}: The array weight should be equal to the average 
#   of the channel (and signal/probe/probeset) weights. Hence, for 
#   single-channel arrays, the array weight should be identical to the 
#   channel weight.
# }
#
#
# \section{Channels and channel weights}{
#   A \emph{channel weight} is a weight assigned to a channel, that is,
#   to the set of signals constituting a channel.
#
#   By definition, a channel weight can never apply across/between arrays.
#
#   \emph{Example}: In two-color data, two channel weights can exist.
#
#   \emph{Example}: In Affymetrix (single-channel) data, only one channel
#   weights can exists and is therefore identical to an array weight.
#
#   \emph{Constraints}: The channel weight should be equal to the average 
#   of all signal/probe/probeset weights in the channel.
# }
#
#
# \section{Combining signal weights into spot weights}{
#   Spot weights can be generated from signal weights. 
#   
#   For a given spot, the spot weight is calculated as the 
#   \emph{arithmetical mean} of the signal weights.
# }
#
#
# \section{Combining signal weights into data-point weights}{
#   Data-point weights can be generated from signal weights. 
#   
#   For a given data point, the data point weight is calculated as the 
#   \emph{arithmetical mean} of the signal weights.
# }
#
#
# \section{Combining spot weights into data-point weights}{
#   Data-point weights can be generated from spot weights. 
#   
#   For a given data point, the data point weight is calculated as the 
#   \emph{arithmetical mean} of the spot weights.
# }
#
# \section{Restrictions}{
#   Note, currently weights are only supported by the @see "RGData" class.
#   The plan is to make this class the "main" class.
#
#   Currently, it is only methods that explicitly say they support
#   weights which use weights. For all other methods, weights are ignored.
# }
#
# @author
#*/#########################################################################


############################################################################
# HISTORY:
# 2005-02-02
# o Clean up.
# 2004-12-19
# o Created.
############################################################################
