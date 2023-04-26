#!/bin/bash

if [[ "$#" -ne 2 ]]; then
    echo ">Usage: mkgmx_fep_decomp [window#-1] [-ligand|-complex]"
    exit 1
fi

if [[ "$2" == "-ligand" ]]; then 
    opt="7 8 9 10 11\n"
elif [[ "$2" == "-complex" ]]; then 
    opt="9 10 11 12 14\n"
else
    echo ">Did not specifiy ligand or complex"
    exit 1
fi

rm -r output
mkdir output

# for jj in $(seq 25 2 47); do for ii in 0 1; do echo -e "9 11 14 $jj $(($jj+1))\n" | gmx energy -f win$ii.edr -o win$ii.xvg >& /dev/null; done;

for ii in $(seq 0 $1)
do
    jj=$((ii+1))
    mkdir output/$ii
    echo -e $opt | gmx energy -f raw/Lambda_${ii}/reswin$ii.edr -o output/$ii/win$ii.xvg >& /dev/null
    echo -e $opt | gmx energy -f raw/Lambda_${ii}/reswin$jj.edr -o output/$ii/win$jj.xvg >& /dev/null
    # cat output/$ii/win${ii}.xvg | grep -v "^@" | grep -v "^#" | awk '{print $8,$3+$6+$7,$2+$4+$5}' >> output/$ii/tmp$ii
    # cat output/$ii/win${ii}.xvg | grep -v "^@" | grep -v "^#" | awk '{print $4,$3,$2}' >> output/$ii/tmp$ii
    cat output/$ii/win${ii}.xvg | grep -v "^@" | grep -v "^#" | awk '{printf "%.3f %.3f %.3f\n", $6,$4+$5,$2+$3}' >> output/$ii/tmp$ii
    cat output/$ii/win${jj}.xvg | grep -v "^@" | grep -v "^#" | awk '{printf "%.3f %.3f %.3f\n", $6,$4+$5,$2+$3}' >> output/$ii/tmp$jj
    paste output/$ii/tmp$ii output/$ii/tmp$jj > output/raw$ii
    # awk '{print ($6 - $1), ($7 + $8 - $2 - $3), ($9 + $10 - $4 - $5)}' < output/raw0 > output/data0
    awk '{print ($4 - $1), ($5 - $2), ($6 - $3)}' < output/raw$ii > output/data$ii

    awk '{esum_elec += (exp($2 * 0.239 / -0.001987 / 298)); \
        esum_vdw += (exp($3 * 0.239 / -0.001987 / 298)); \
        esum_total += (exp($1 * 0.239 / -0.001987 / 298)); \
        sum_elec += $2 * 0.239; sum2_elec += $2*$2 * 0.057; \
        sum_vdw += $3 * 0.239; sum2_vdw += $3*$3 * 0.057; \
        sum_couple += $2*$3 * 0.057; n++;} END \
        {printf "%3d % 6.2f % 10.2f % 6.2f % 6.2f % 6.2f % 6.2f % 10.2f % 6.2f % 6.2f % 6.2f % 6.2f % 10.2f\n", \
            n, \
            (-0.001987*298*log(esum_total/n)), \
            (-0.001987*298*log(esum_elec/n)), \
            (sum_elec/n) - (sum2_elec/n - sum_elec*sum_elec/n/n)/1.24 - (sum_couple/n - sum_elec*sum_vdw/n/n)/1.24, \
            (sum_elec/n), \
            -(sum2_elec/n - sum_elec*sum_elec/n/n)/1.24, \
            -(sum_couple/n - sum_elec*sum_vdw/n/n)/1.24, \
            (-0.001987*298*log(esum_vdw/n)), \
            (sum_vdw/n) - (sum2_vdw/n - sum_vdw*sum_vdw/n/n)/1.24 - (sum_couple/n - sum_elec*sum_vdw/n/n)/1.24, \
            (sum_vdw/n), \
            -(sum2_vdw/n - sum_vdw*sum_vdw/n/n)/1.24, \
            -(sum_couple/n - sum_elec*sum_vdw/n/n)/1.24, \
            -0.001987*298*log(esum_total/n) + 0.001987*298*log(esum_elec/n) + 0.001987*298*log(esum_vdw/n) \
            }' < output/data$ii >> output/decomp
    awk '{sum1 += (exp($1 * 0.239 / -0.001987 / 298)); \
        sum2 += (exp(2 * $1 * 0.239 / -0.001987 / 298)); \
        sum_total += $1 * 0.239; n++;} END \
        {printf "%3d % 6.2f %5.2f % 6.2f %5.2f %7.0f %3.2e %3.2e %3.2e\n", \
            n, \
            (-0.001987*298*log(sum1/n)), \
            sqrt(((n*sum2/sum1/sum1-1)/0.3844)), \
            (sum_total/n), \
            (sum_total/n + 0.001987*298*log(sum1/n)), \
            (exp((sum_total/n - (-0.001987*298*log(sum1/n)))/298/0.001987)), \
            (sum1/n), \
            (sum1*sum1/n/n), \
            (sum2/n)}' < output/data$ii >> output/error
done

# echo -n ">>>>FEP gives you: "
# tail -1 $1 | awk '{print $19}'
awk '{sum_dg += $2; \
    sum_delec += $3; \
    sum_dvdw += $8; \
    sum_couple += $13; n++;} END \
    {printf ">Decomp gives you: %.2f\nin which \nelec : % 7.2f\nvdw  : % 7.2f\ncoupl: % 7.2f\n", sum_dg, sum_delec, sum_dvdw, sum_couple}' < output/decomp
