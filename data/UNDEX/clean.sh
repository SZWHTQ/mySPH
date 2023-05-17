#!/usr/bin/env zsh


if [ ! -d "./output" ]
then
    mkdir -p ./output
fi

if [ ! -d "./vtk" ]
then
    mkdir -p ./vtk
fi

# Clean up the data directory
if [ "`exa ./output`" != "" ]
then
    rm ./output/*.dat
fi

if [ "`exa ./vtk`" != "" ]
then
    rm ./vtk/*.vtk
fi

if [ -f "./ini.vtk" ]
then
    rm ./ini.vtk
fi


if [ "`exa -a ./output`" != "" ]
then
   echo "Directory: './output' still not empty"
   exa -a ./output/
fi