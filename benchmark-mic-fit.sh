#! /bin/bash

sed -i 's/int nTracks = 20000/int nTracks = 1000000/g' Config.cc 

make -j 8

dir=/data/nfsmic/${USER}/tmp
micdir=/nfsmic/${USER}/tmp

mkdir -p ${dir}
./mkFit/mkFit --write --file-name simtracks_10x1M.bin
mv simtracks_10x1M.bin ${dir}/

for nth in 1 3 7 21 42 63 84 105 126 147 168 189 210
do
echo nth=${nth}
ssh mic0 ./mkFit-mic --read --file-name ${micdir}/simtracks_10x1M.bin --fit-std-only --num-thr ${nth} >& log_mic_10x1M_FIT_NVU16int_NTH${nth}.txt
done

sed -i 's/# USE_INTRINSICS := -DMPT_SIZE=1/USE_INTRINSICS := -DMPT_SIZE=XX/g' Makefile.config
for nvu in 1 2 4 8 16
do
sed -i "s/MPT_SIZE=XX/MPT_SIZE=${nvu}/g" Makefile.config
make clean
make -j 8
echo nvu=${nvu}
ssh mic0 ./mkFit-mic --read --file-name ${micdir}/simtracks_10x1M.bin --fit-std-only --num-thr 1 >& log_mic_10x1M_FIT_NVU${nvu}_NTH1.txt
sed -i "s/MPT_SIZE=${nvu}/MPT_SIZE=XX/g" Makefile.config
done
sed -i 's/USE_INTRINSICS := -DMPT_SIZE=XX/# USE_INTRINSICS := -DMPT_SIZE=1/g' Makefile.config

sed -i 's/int nTracks = 1000000/int nTracks = 20000/g' Config.cc 

make clean
make -j 8
