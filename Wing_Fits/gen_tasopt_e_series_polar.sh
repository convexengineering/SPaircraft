#airfoil=$1
POLARFILE=blade.e$1.Re$2k.M$3.pol

if [ -f $POLARFILE ] ; then
    echo "yes"
    rm $POLARFILE
fi

xfoil << EOF
load blade.e$1
pane blade.e$1
oper
v $2e3
M 0.1
a 0
M 0.15
a 0
M 0.2
a 0
M 0.25
a 0
M 0.3
a 0
M 0.35
a 0
M 0.4
a 0
M 0.45
a 0
M 0.5
a 0
M 0.55
a 0
M 0.6
a 0
M 0.65
a 0
M 0.7
a 0
M 0.75
a 0
M 0.8
a 0
M 0.85
a 0
M 0.9
a 0
M $3
pacc
$POLARFILE

iter 100
cseq .3 .7 .05

quit
EOF
