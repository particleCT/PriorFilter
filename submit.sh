#!/bin/sh                                                                                                                                                                                                                                   #PBS -N binary                                                                                                                                                                                                                              #PBS -e output.err                                                                                                                                                                                                                          #PBS -o output.out 

source /home/cacfek/.bashrc
cd /ion/home/hect/Filter/PriorFilter

/ion/home/hect/Filter/PriorFilter/bin/PriorFilter /ion/home/hect/Prior/CIRSPHP1/inf/Prior_FBP.root /ion/home/hect/DataConversion/Bin2Root/Data/projection_307.0.root &
/ion/home/hect/Filter/PriorFilter/bin/PriorFilter /ion/home/hect/Prior/CIRSPHP1/inf/Prior_FBP.root /ion/home/hect/DataConversion/Bin2Root/Data/projection_308.5.root &
/ion/home/hect/Filter/PriorFilter/bin/PriorFilter /ion/home/hect/Prior/CIRSPHP1/inf/Prior_FBP.root /ion/home/hect/DataConversion/Bin2Root/Data/projection_126.5.root &
/ion/home/hect/Filter/PriorFilter/bin/PriorFilter /ion/home/hect/Prior/CIRSPHP1/inf/Prior_FBP.root /ion/home/hect/DataConversion/Bin2Root/Data/projection_273.5.root &
/ion/home/hect/Filter/PriorFilter/bin/PriorFilter /ion/home/hect/Prior/CIRSPHP1/inf/Prior_FBP.root /ion/home/hect/DataConversion/Bin2Root/Data/projection_70.5.root &
/ion/home/hect/Filter/PriorFilter/bin/PriorFilter /ion/home/hect/Prior/CIRSPHP1/inf/Prior_FBP.root /ion/home/hect/DataConversion/Bin2Root/Data/projection_214.0.root &
/ion/home/hect/Filter/PriorFilter/bin/PriorFilter /ion/home/hect/Prior/CIRSPHP1/inf/Prior_FBP.root /ion/home/hect/DataConversion/Bin2Root/Data/projection_47.0.root &
/ion/home/hect/Filter/PriorFilter/bin/PriorFilter /ion/home/hect/Prior/CIRSPHP1/inf/Prior_FBP.root /ion/home/hect/DataConversion/Bin2Root/Data/projection_60.5.root &
wait
