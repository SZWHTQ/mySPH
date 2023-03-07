#!/usr/bin/env zsh

# Clean and complete this project

if [ ! -d "./data/demo/shockTube/output" ]
then
    mkdir -p ./data/demo/shockTube/output
fi

if [ ! -d "./data/demo/shockTube/vtk" ]
then
    mkdir -p ./data/demo/shockTube/vtk
fi

if [ ! -d "./data/demo/shearCavity/output" ]
then
    mkdir -p ./data/demo/shearCavity/output
fi

if [ ! -d "./data/demo/shearCavity/vtk" ]
then
    mkdir -p ./data/demo/shearCavity/vtk
fi

if [ ! -d "./data/tntBar/output" ]
then
    mkdir -p ./data/tntBar/output
fi

if [ ! -d "./data/tntBar/vtk" ]
then
    mkdir -p ./data/tntBar/vtk
fi

if [ ! -d "./data/tntBar/pic" ]
then
    mkdir -p ./data/tntBar/pic
fi

if [ ! -d "./data/tntCylinder/output" ]
then
    mkdir -p ./data/tntCylinder/output
fi

if [ ! -d "./data/tntCylinder/vtk" ]
then
    mkdir -p ./data/tntCylinder/vtk
fi

if [ ! -d "./data/undexCylinder/output" ]
then
    mkdir -p ./data/undexCylinder/output
fi

if [ ! -d "./data/undexCylinder/vtk" ]
then
    mkdir -p ./data/undexCylinder/vtk
fi

if [ ! -d "./data/undexChamber/output" ]
then
    mkdir -p ./data/undexChamber/output
fi

if [ ! -d "./data/undexChamber/vtk" ]
then
    mkdir -p ./data/undexChamber/vtk
fi


if [ "`exa ./data/demo/shearCavity/output`" != "" ]
then
    rm ./data/demo/shearCavity/output/*.dat
fi

if [ "`exa ./data/demo/shearCavity/vtk`" != "" ]
then
    rm ./data/demo/shearCavity/vtk/*.vtk
fi

if [ "`exa ./data/demo/shockTube/output`" != "" ]
then
    rm ./data/demo/shockTube/output/*.dat
fi

if [ "`exa ./data/demo/shockTube/vtk`" != "" ]
then
    rm ./data/demo/shockTube/vtk/*.vtk
fi

if [ "`exa ./data/tntBar/output`" != "" ]
then
    rm ./data/tntBar/output/*.dat
fi

if [ "`exa ./data/tntBar/vtk`" != "" ]
then
    rm ./data/tntBar/vtk/*.vtk
fi

if [ "`exa ./data/tntBar/pic`" != "" ]
then
    rm ./data/tntBar/pic/*.png
fi

if [ "`exa ./data/tntCylinder/output`" != "" ]
then
    rm ./data/tntCylinder/output/*.dat
fi

if [ "`exa ./data/tntCylinder/vtk`" != "" ]
then
    rm ./data/tntCylinder/vtk/*.vtk
fi

if [ "`exa ./data/undexCylinder/output`" != "" ]
then
    rm ./data/undexCylinder/output/*.dat
fi

if [ "`exa ./data/undexCylinder/vtk`" != "" ]
then
    rm ./data/undexCylinder/vtk/*.vtk
fi

if [ "`exa ./data/undexChamber/output`" != "" ]
then
    rm ./data/undexChamber/output/*.dat
fi

if [ "`exa ./data/undexChamber/vtk`" != "" ]
then
    rm ./data/undexChamber/vtk/*.vtk
fi

if [ -f "./data/demo/shearCavity/ini.vtk" ] || [ -f "./data/demo/shockTube/ini.vtk" ]
then
    rm ./data/demo/*/ini.vtk
fi

if [ -f "./data/tntBar/ini.vtk" ] || [ -f "./data/tntCylinder/ini.vtk" ] || [ -f "./data/undexCylinder/ini.vtk" ] || [ -f "./data/undexChamber/ini.vtk" ]
then
    rm ./data/*/ini.vtk
fi

if [ "`find ./ -name '._*'`" != "" ]
then
    find ./ -name "._*" | xargs rm
fi

if [ "`find ./ -name '.DS_*'`" != "" ]
then
    rm `find ./ -name ".DS_*"`
fi


if [ "`exa -a ./data/demo/shearCavity/output`" != "" ]
then
   echo "Directory: './data/demo/shearCavity/output' still not empty"
   exa -a ./data/demo/shearCavity/output/
fi

if [ "`exa -a ./data/demo/shearCavity/vtk`" != "" ]
then
   echo "Directory: './data/demo/shearCavity/vtk' still not empty"
   exa -a ./data/demo/shearCavity/vtk/
fi

if [ "`exa -a ./data/demo/shockTube/output`" != "" ]
then
   echo "Directory: './data/demo/shockTube/output' still not empty"
   exa -a ./data/demo/shockTube/output/
fi

if [ "`exa -a ./data/demo/shockTube/vtk`" != "" ]
then
   echo "Directory: './data/demo/shockTube/vtk' still not empty"
   exa -a ./data/demo/shockTube/vtk/
fi

if [ "`exa -a ./data/tntBar/output`" != "" ]
then
   echo "Directory: './data/tntBar/output' still not empty"
   exa -a ./data/tntBar/output/
fi

if [ "`exa -a ./data/tntBar/vtk`" != "" ]
then
   echo "Directory: './data/tntBar/vtk' still not empty"
   exa -a ./data/tntBar/vtk/
fi

if [ "`exa -a ./data/tntCylinder/output`" != "" ]
then
   echo "Directory: './data/tntCylinder/output' still not empty"
   exa -a ./data/tntCylinder/output/
fi

if [ "`exa -a ./data/tntCylinder/vtk`" != "" ]
then
   echo "Directory: './data/tntCylinder/vtk' still not empty"
   exa -a ./data/tntCylinder/vtk/
fi

if [ "`exa -a ./data/undexCylinder/output`" != "" ]
then
   echo "Directory: './data/undexCylinder/output' still not empty"
   exa -a ./data/undexCylinder/output/
fi

if [ "`exa -a ./data/undexCylinder/vtk`" != "" ]
then
   echo "Directory: './data/undexCylinder/vtk' still not empty"
   exa -a ./data/undexCylinder/vtk/
fi

if [ "`exa -a ./data/undexChamber/output`" != "" ]
then
   echo "Directory: './data/undexChamber/output' still not empty"
   exa -a ./data/undexChamber/output/
fi

if [ "`exa -a ./data/undexChamber/vtk`" != "" ]
then
   echo "Directory: './data/undexChamber/vtk' still not empty"
   exa -a ./data/undexChamber/vtk/
fi

./data/damBreak/clean.sh

if [ $# != 0 ]
then
    if [ $1 = "--all" -o $1 = "-all" ]
    then
        if [ -d "./build" ] && [ "`exa ./build`" != "" ]
        then
           rm -rf ./build/*
        fi

        if [ -d "./build" ] && [ "`exa -a ./build`" != "" ]
        then
            rm -rf ./build/.*
        fi

        if [ -d "./_build" ]
        then
           rm -rf ./_build
        fi

        if [ "`exa ./module`" != "" ]
        then
            rm -rf module/*.mod
        fi
        
        if [ -d "./dependencies" ]
        then
            rm -rf ./dependencies
        fi

    fi
fi