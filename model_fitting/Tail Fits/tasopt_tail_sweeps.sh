
thick="100 120 140"
Re="500 1000 1500 2000 2500 3000 3500 4000 4500 5000 5500 6000 6500 7000 7500 8000 8500 9000"
M="0.4 0.6 0.8"

for m in $M
do
    for r in $Re
    do
        for t in $thick
        do
            ./gen_tasopt_tail_polar.sh $t $r $m
        done
    done
done
