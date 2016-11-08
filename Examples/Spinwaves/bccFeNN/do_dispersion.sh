#!/bin/bash
echo "Performing time-series"
./timeseries kvecinfo.dat 1e12 1e10 5e15
filename=bccFe_SW_Gamma_H_N_Gamma.dat
if [ -f $filename ];then
    echo "$filename removed, starting again."
    rm $filename
fi
nkp=`grep NumberKPoints kvecinfo.dat | awk '{print $2-1}'`
echo "Number of k-points"
for i in `seq 0 1 $nkp`;do
    #get the point
    kx=`grep "lookupkvector $i =" kvecinfo.dat | awk '{print $4}'`
    ky=`grep "lookupkvector $i =" kvecinfo.dat | awk '{print $5}'`
    kz=`grep "lookupkvector $i =" kvecinfo.dat | awk '{print $6}'`
    echo "k-point $i, kx = $kx , ky = $ky , kz = $kz"
    if [ -f kx${kx}ky${ky}kz${kz}.dat ];then
        cat kx${kx}ky${ky}kz${kz}.dat | awk -v k=$i -v kx=$kx -v ky=$ky -v kz=$kz '{print k,kx,ky,kz,$1,$2,$3,$4,$5}' >> $filename
        echo "" >> $filename
        echo "" >> $filename
    else
        echo "ERROR: Could not find kx${kx}ky${ky}kz${kz}.dat"
    fi
done
