# setupPlots.sh
#
# David Adams
# May 2019
#
# Example configuration file for the drawNoise command in dunenoise.

SDET=iceberg3
PLOTS="tai-nsgrms:300,tai-nsgrms50:200,cnr-nsgrms:300,cnr-nsgrms50:200"
CRVS="zcGood-uvGood"
DIRPAT="ib%REC%NoiseVsChan/iceberg%RUN%/chmet_%VAR%_all_run%RUN%.tpad"
GRPPAT="ib%REC%NoiseVsChan/ibgroup%RUN%/chmet_%VAR%_all_run%RUN1%.tpad"
OUTPRE="plots-noisesum/"
