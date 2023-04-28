#!/usr/bin/env zsh

# Clean and complete this project

if [ $# != 0 ]
then
    if [ $1 = "--all" -o $1 = "-all" ]
    then
        cd data
        ./clean.sh
        cd ..
    fi
fi

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