#! /bin/bash
runtag(){
#tags='full sigall sigzpt100 sigpt100'
tags='sigzpt50_log sigpt100_log'
#tags=zveto


for itag in $tags
do
    echo "command: $1 $itag"
    python $1 $itag >$1_$itag.log &
done
}
if [ -z "$1" ]
then
    echo -e "\n$0 script1 [script2]  [script3] ...\n"
    exit
fi
for whatscript in $*
do
    runtag $whatscript
done

