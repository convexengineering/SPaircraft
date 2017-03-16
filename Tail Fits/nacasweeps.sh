
NACA="0005 0008 0009 0010 0015 0020"
Re="500 1000 1500 2000 2500 3000 3500 4000 4500 5000 5500 6000 6500 7000 7500 8000 8500 9000"

for r in $Re
do
    for n in $NACA
    do
        ./gen_tail_polar.sh $n $r
    done
done
