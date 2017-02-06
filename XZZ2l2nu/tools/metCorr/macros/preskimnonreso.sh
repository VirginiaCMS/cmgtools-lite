#!/bin/sh



# compile
g++ preskimnonreso.cc -o preskimnonreso.exe `root-config --cflags` `root-config --libs`


./preskimnonreso.exe /data2/XZZ2/80X_20170124_light_Skim/muonegtree_light_skim_38.root /home/heli/XZZ/80X_20170202_light_Skim/muonegtree_light_skim_38_skim.root 
