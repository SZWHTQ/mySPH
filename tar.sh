#!/usr/bin/env zsh
./clean.sh --all
name="mySPH_`date "+%Y%m%d-%H%M"`.tgz"
cd ..
tar --exclude doc/books/光滑粒子流体动力学：一种无网格粒子法.pdf \
    --use-compress-program=pigz \
    -pcvf ~/git/TimeMachine/$name ./mySPH
chmod 755 ~/git/TimeMachine/$name
cd ./mySPH