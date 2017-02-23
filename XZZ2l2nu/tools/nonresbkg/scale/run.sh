#! /bin/bash

cuts='0 10 20 30 40 50 60 70 80 90 100 110'


echo -e "ZPT cut\t\t\t   electron\t   \t       muon">log
for acut in $cuts
do
    python histmaker.py el $acut &
    python histmaker.py mu $acut &
    wait
    esc=$(python getscale elconfig)
    msc=$(python getscale muconfig)
    echo -e "$acut GeV   \t$esc\t   $msc" >>log
done