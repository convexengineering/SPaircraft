
thick="090 100 110 120 130 140 145"
#Re="500 1000 1500 2000 2500 3000 3500 4000 4500 5000 5500 6000 6500 7000 7500 8000 8500 9000"
Re="10000 15000 20000 25000 30000"
M="0.4 0.5 0.6 0.7 0.8 0.9"

for m in $M
do
    for r in $Re
    do
        for t in $thick
        do
            ./gen_tasopt_c_series_polar.sh $t $r $m
        done
    done
done
