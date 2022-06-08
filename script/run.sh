#!/bin/bash

echo ==========================================================
echo = Executing channel2d on all files inside ../data/channels
echo ==========================================================

for file in ../data/channels/*.chn
do
    output1=`echo $file | sed -e 's/channels/spcurves/'`
    output2=`echo $file | sed -e 's/channels/lineprogs/'`
    output1=`echo $output1 | sed -e 's/\.chn/\.spl/'`
    output2=`echo $output2 | sed -e 's/\.chn/\.pl/'`
  
    echo executing ../bin/Channel2D-App $file $output1 $output2
    echo `../bin/Channel2D-App $file $output1 $output2`
done
