param l {i in 0..1};
param u {i in 0..1};
param xini {i in 0..1};

var x {i in 0..1} >= l[i], <= u[i], := xini[i];

minimize Obj:
    x[1];

subject to R1:
    x[0]^2 + 1 - x[1] <= 0;

subject to R2:
    2 - x[0] - x[1] <= 0;

data;

param xini :=
    0   0
    1   0
;

param l :=
    0   -10
    1   -1.0e+20
;

param u :=
    0    10
    1    1.0e+20
;

