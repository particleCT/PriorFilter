#!/bin/sh                                                                                                                                                                                                                                   #PBS -N binary                                                                                                                                                                                                                              #PBS -e output.err                                                                                                                                                                                                                          #PBS -o output.out 

source /home/cacfek/.bashrc
cd /ion/home/hect/Filter/PriorFilter

$COMMANDLINE
wait
