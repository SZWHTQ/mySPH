#!/usr/bin/env zsh
type0_files=(output/Type_0*.dat)
type5_files=(output/Type_5*.dat)

for i in {0..$#type5_files}; do
  cat $type5_files[$i] <(tail -n +2 $type0_files[$i]) > "output/Reconstruct_$i.dat"
done
