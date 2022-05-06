OPENQASM 2.0;
include "qelib1.inc";
qreg a[2];
qreg b[1];
id a[0];
id a[1];
id b[0];
rx(1.9948058) a[0];
ry(0.17145549) a[0];
rz(1.194034) a[0];
rx(0.42890442) a[1];
ry(1.5978195) a[1];
rz(0.45533227) a[1];
barrier a[0],a[1],b[0];
id b[0];
h b[0];
id a[0];
h a[0];
id a[1];
sdg a[1];
cx b[0],a[1];
cx b[0],a[0];
h a[0];
s a[1];
barrier a[0],a[1],b[0];
id b[0];
id a[0];
s a[0];
id a[1];
sdg a[1];
h a[1];
cx a[0],a[1];
h a[0];
s a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
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
