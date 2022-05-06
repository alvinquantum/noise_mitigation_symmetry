OPENQASM 2.0;
include "qelib1.inc";
qreg a[2];
qreg b[1];
id a[0];
id a[1];
id b[0];
rx(0.36693431) a[0];
ry(1.9490507) a[0];
rz(0.60787599) a[0];
rx(0.8773267) a[1];
ry(1.9828644) a[1];
rz(1.1787228) a[1];
barrier a[0],a[1],b[0];
id b[0];
h b[0];
id a[0];
id a[1];
cx b[0],a[1];
cx b[0],a[0];
barrier a[0],a[1],b[0];
id b[0];
id a[0];
y a[0];
h a[0];
id a[1];
y a[1];
h a[1];
s a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
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
