#!/bin/sh



# compile
g++ preskimnonreso.cc -o preskimnonreso.exe `root-config --cflags` `root-config --libs`


./preskimnonreso.exe /data2/yanchu/muoneg36p4/muonegtree.root  /home/heli/XZZ/80X_20161029_light/muonegtree_light.root 
