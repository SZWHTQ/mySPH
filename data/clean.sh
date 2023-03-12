#!/usr/bin/env zsh

for file in ./*
do
    if [ -d "$file" ]
    then
        cd "$file"
        if [ -f "./clean.sh" ]
        then
            ./clean.sh
        else
            for file in ./*
            do
                if [ -d "$file" ]
                then
                    cd "$file"
                    if [ -f "./clean.sh" ]
                    then
                        ./clean.sh
                    fi
                    cd ..
                    echo "Cleaned $file"
                fi
            done
        fi
        cd ..
        echo "Cleaned $file"
    fi
done