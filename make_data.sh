#!/bin/zsh

HOME=`pwd`

# EBL
cd $HOME/ebl
python convert.py
mv ebl_*.txt $HOME/data

# LOSSES
cd $HOME/losses
python convert.py
mv losses_*.txt $HOME/data

# XSECS
cd $HOME/xsecs
python convert.py
mv xsecs_*.txt $HOME/data
