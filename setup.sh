# setup.sh
#
# David Adams
# January 2020

if [ -n "$DUNENOISE_DIR" ]; then
  echo It appears dunenoise is already set up.
  echo DUNENOISE_DIR = $DUNENOISE_DIR

elif [ -z "$DUNEPROC_DIR" ]; then
  echo First set up duneproc

else

echo Setting up dunenoise
export DUNENOISE_DIR=$(dirname $(readlink -f $BASH_SOURCE))
PATH=$DUNENOISE_DIR/bin:$PATH
TMPPATH=.:./job:
if [ ${FHICL_FILE_PATH:0:8} = $TMPPATH ]; then FHICL_FILE_PATH=${FHICL_FILE_PATH:8}; fi
FHICL_FILE_PATH=${TMPPATH}$DUNENOISE_DIR/fcl:$FHICL_FILE_PATH
unset TMPPATH

fi
