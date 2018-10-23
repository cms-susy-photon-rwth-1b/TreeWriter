#####!/bin/bash

#export SCRAM_ARCH=slc6_amd64_gcc530
export SCRAM_ARCH=slc7_amd64_gcc530
export VO_CMS_SW_DIR=/cvmfs/cms.cern.ch
source $VO_CMS_SW_DIR/cmsset_default.sh



. /cvmfs/cms.cern.ch/crab3/crab.sh

for jobDir in crab_*
   do
      echo "killing " $jobDir
      crab kill $jobDir
      echo "remove folder"
      rm -r $jobDir
      echo ""
   done
