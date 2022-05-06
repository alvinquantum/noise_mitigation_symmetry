OPENQASM 2.0;
include "qelib1.inc";
qreg a[2];
qreg b[3];
id a[0];
id a[1];
id b[0];
id b[1];
id b[2];
rx(0.39519852) a[0];
ry(0.7010435) a[0];
rz(1.0617151) a[0];
rx(1.2413007) a[1];
ry(0.36982502) a[1];
rz(1.2288855) a[1];
barrier a[0],a[1],b[0],b[1],b[2];
id a[0];
sdg a[0];
id a[1];
id b[0];
h b[0];
id b[1];
h b[1];
cx b[1],a[1];
h a[1];
cx b[0],a[1];
h a[1];
id b[2];
h b[2];
cx b[2],a[0];
s a[0];
sdg a[0];
cx b[1],a[0];
s a[0];
sdg a[0];
cx b[0],a[0];
s a[0];
barrier a[0],a[1],b[0],b[1],b[2];
id a[0];
h a[0];
s a[0];
h a[0];
id a[1];
x a[1];
s a[1];
y a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
sdg a[1];
h a[1];
id b[0];
id b[1];
id b[2];
barrier a[0],a[1],b[0],b[1],b[2];
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
id b[2];
cx b[2],a[0];
h b[2];
barrier a[0],a[1],b[0],b[1],b[2];
id a[0];
id a[1];
id b[0];
id b[1];
id b[2];
