OPENQASM 2.0;
include "qelib1.inc";
qreg a[2];
qreg b[2];
id a[0];
id a[1];
id b[0];
id b[1];
rx(1.1319868) a[0];
ry(1.5503016) a[0];
rz(1.4881774) a[0];
rx(0.28095804) a[1];
ry(1.0709049) a[1];
rz(1.0679557) a[1];
barrier a[0],a[1],b[0],b[1];
id a[0];
id a[1];
h a[1];
id b[0];
h b[0];
id b[1];
h b[1];
cx b[1],a[1];
h a[1];
sdg a[1];
cx b[0],a[1];
s a[1];
cx b[1],a[0];
cx b[0],a[0];
barrier a[0],a[1],b[0],b[1];
id a[0];
y a[0];
h a[0];
id a[1];
y a[1];
s a[1];
h a[1];
s a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
sdg a[1];
h a[1];
id b[0];
id b[1];
barrier a[0],a[1],b[0],b[1];
id a[0];
id a[1];
id b[0];
cx b[0],a[1];
h a[1];
h b[0];
id b[1];
cx b[1],a[1];
h a[1];
h b[1];
barrier a[0],a[1],b[0],b[1];
id a[0];
id a[1];
id b[0];
id b[1];
