OPENQASM 2.0;
include "qelib1.inc";
qreg a[2];
qreg b[1];
barrier a[0],a[1],b[0];
id a[0];
id a[1];
id b[0];
h b[0];
sdg a[1];
cx b[0],a[1];
s a[1];
cx b[0],a[0];
barrier a[0],a[1],b[0];
id a[0];
id a[1];
id b[0];
z a[0];
h a[0];
s a[0];
h a[1];
cx a[0],a[1];
h a[0];
s a[0];
h a[1];
s a[1];
cx a[0],a[1];
barrier a[0],a[1],b[0];
id a[0];
id a[1];
id b[0];
cx b[0],a[1];
h b[0];
barrier a[0],a[1],b[0];
id a[0];
id a[1];
id b[0];
