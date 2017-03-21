#airfoil=$1
POLARFILE=blade.t$1.cl0.Re$2k.M$3.pol

if [ -f $POLARFILE ] ; then
    echo "yes"
    rm $POLARFILE
fi

xfoil << EOF
load blade.t$1
pane blade.t$1
oper
v $2e3
M $3
pacc
$POLARFILE

iter 400
cl 0.0

quit
EOF
