OPENQASM 2.0;
include "qelib1.inc";
qreg a[2];
qreg b[1];
id a[0];
id a[1];
id b[0];
rx(0.58510236) a[0];
ry(0.19819513) a[0];
rz(0.74849156) a[0];
rx(1.8659579) a[1];
ry(1.6757217) a[1];
rz(1.5964202) a[1];
barrier a[0],a[1],b[0];
id b[0];
h b[0];
id a[0];
id a[1];
sdg a[1];
cx b[0],a[1];
s a[1];
barrier a[0],a[1],b[0];
id b[0];
id a[0];
x a[0];
s a[0];
id a[1];
z a[1];
cx a[0],a[1];
h a[0];
s a[0];
x a[0];
sdg a[0];
h a[0];
sdg a[1];
h a[1];
z a[1];
h a[1];
s a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
h a[1];
s a[1];
barrier a[0],a[1],b[0];
id b[0];
id a[0];
id a[1];
cx b[0],a[1];
h b[0];
barrier a[0],a[1],b[0];
id b[0];
id a[0];
id a[1];
