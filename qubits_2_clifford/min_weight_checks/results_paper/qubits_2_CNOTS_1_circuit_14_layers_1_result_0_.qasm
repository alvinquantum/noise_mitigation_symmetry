OPENQASM 2.0;
include "qelib1.inc";
qreg a[2];
qreg b[1];
id a[0];
id a[1];
id b[0];
rx(0.29993596) a[0];
ry(0.31741116) a[0];
rz(1.2565144) a[0];
rx(0.34578284) a[1];
ry(1.4654707) a[1];
rz(0.87609676) a[1];
barrier a[0],a[1],b[0];
id b[0];
h b[0];
id a[0];
h a[0];
id a[1];
cx b[0],a[1];
cx b[0],a[0];
h a[0];
barrier a[0],a[1],b[0];
id b[0];
id a[0];
id a[1];
x a[1];
h a[1];
s a[1];
cx a[0],a[1];
h a[0];
s a[0];
sdg a[1];
h a[1];
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
