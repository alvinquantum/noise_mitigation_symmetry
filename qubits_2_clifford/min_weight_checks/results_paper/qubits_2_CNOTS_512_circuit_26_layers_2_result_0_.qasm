OPENQASM 2.0;
include "qelib1.inc";
qreg a[2];
qreg b[2];
id a[0];
id a[1];
id b[0];
id b[1];
rx(0.75442374) a[0];
ry(0.30574872) a[0];
rz(1.4930766) a[0];
rx(0.60040633) a[1];
ry(0.58837693) a[1];
rz(1.656881) a[1];
barrier a[0],a[1],b[0],b[1];
id a[0];
sdg a[0];
id a[1];
id b[0];
h b[0];
cx b[0],a[1];
id b[1];
h b[1];
cx b[1],a[0];
s a[0];
cx b[0],a[0];
barrier a[0],a[1],b[0],b[1];
id a[0];
x a[0];
s a[0];
z a[0];
s a[0];
h a[0];
s a[0];
id a[1];
y a[1];
h a[1];
s a[1];
y a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
x a[0];
s a[0];
z a[1];
h a[1];
s a[1];
cx a[0],a[1];
h a[0];
s a[0];
h a[1];
s a[1];
cx a[0],a[1];
y a[0];
h a[0];
x a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
h a[0];
s a[0];
h a[1];
s a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
z a[1];
h a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
h a[0];
s a[0];
y a[0];
h a[1];
s a[1];
cx a[0],a[1];
h a[0];
s a[0];
h a[1];
s a[1];
cx a[0],a[1];
x a[0];
s a[0];
z a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
h a[1];
s a[1];
cx a[0],a[1];
y a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
h a[0];
s a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
y a[0];
s a[0];
h a[0];
s a[0];
z a[1];
h a[1];
s a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
z a[0];
h a[0];
s a[0];
s a[1];
cx a[0],a[1];
h a[0];
s a[0];
h a[0];
sdg a[1];
h a[1];
x a[1];
s a[1];
h a[1];
s a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
z a[0];
h a[0];
s a[0];
h a[1];
s a[1];
x a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
x a[0];
s a[0];
h a[0];
s a[0];
y a[1];
h a[1];
s a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
h a[1];
s a[1];
cx a[0],a[1];
z a[0];
h a[0];
z a[1];
cx a[0],a[1];
h a[0];
s a[0];
x a[0];
sdg a[0];
h a[0];
z a[1];
s a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
y a[0];
s a[0];
h a[1];
s a[1];
x a[1];
h a[1];
s a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
x a[0];
s a[0];
h a[0];
s a[0];
x a[0];
s a[0];
sdg a[1];
h a[1];
sdg a[1];
h a[1];
z a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
z a[0];
s a[0];
y a[1];
s a[1];
cx a[0],a[1];
h a[0];
s a[0];
z a[0];
sdg a[1];
h a[1];
x a[1];
s a[1];
cx a[0],a[1];
h a[0];
s a[0];
y a[0];
s a[0];
h a[0];
s a[0];
sdg a[1];
h a[1];
s a[1];
cx a[0],a[1];
h a[0];
s a[0];
sdg a[0];
h a[0];
z a[1];
s a[1];
h a[1];
s a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
x a[0];
sdg a[0];
h a[0];
s a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
h a[0];
s a[0];
h a[0];
s a[0];
h a[1];
s a[1];
z a[1];
s a[1];
cx a[0],a[1];
h a[0];
s a[0];
y a[0];
h a[0];
h a[1];
s a[1];
sdg a[1];
h a[1];
cx a[0],a[1];
h a[0];
s a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
z a[0];
s a[0];
h a[1];
s a[1];
z a[1];
sdg a[1];
h a[1];
cx a[0],a[1];
h a[0];
s a[0];
y a[0];
h a[1];
s a[1];
y a[1];
h a[1];
s a[1];
cx a[0],a[1];
h a[0];
s a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
z a[0];
h a[0];
y a[1];
h a[1];
cx a[0],a[1];
h a[0];
s a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
y a[0];
h a[1];
s a[1];
x a[1];
sdg a[1];
h a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
h a[0];
s a[0];
y a[0];
s a[0];
h a[1];
s a[1];
z a[1];
h a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
y a[0];
h a[0];
sdg a[1];
h a[1];
y a[1];
cx a[0],a[1];
h a[0];
s a[0];
z a[0];
s a[0];
h a[0];
s a[0];
h a[1];
s a[1];
h a[1];
s a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
h a[0];
s a[0];
z a[0];
h a[0];
h a[1];
s a[1];
s a[1];
cx a[0],a[1];
y a[0];
h a[0];
h a[1];
s a[1];
h a[1];
s a[1];
cx a[0],a[1];
z a[0];
s a[0];
h a[0];
s a[0];
h a[1];
s a[1];
s a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
sdg a[1];
h a[1];
y a[1];
s a[1];
h a[1];
s a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
x a[0];
h a[0];
s a[0];
h a[1];
s a[1];
z a[1];
s a[1];
h a[1];
s a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
h a[0];
s a[0];
h a[1];
s a[1];
cx a[0],a[1];
y a[0];
h a[0];
s a[0];
z a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
z a[0];
z a[1];
sdg a[1];
h a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
h a[0];
s a[0];
x a[0];
h a[0];
s a[0];
h a[1];
s a[1];
y a[1];
s a[1];
h a[1];
s a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
y a[0];
h a[1];
s a[1];
x a[1];
s a[1];
h a[1];
s a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
h a[0];
s a[0];
h a[1];
s a[1];
cx a[0],a[1];
s a[0];
h a[0];
s a[0];
y a[1];
h a[1];
s a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
y a[0];
sdg a[0];
h a[0];
h a[1];
s a[1];
x a[1];
h a[1];
s a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
h a[0];
s a[0];
z a[0];
s a[0];
z a[1];
h a[1];
cx a[0],a[1];
z a[0];
sdg a[1];
h a[1];
z a[1];
h a[1];
s a[1];
cx a[0],a[1];
h a[1];
s a[1];
y a[1];
sdg a[1];
h a[1];
cx a[0],a[1];
h a[0];
s a[0];
z a[0];
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
cx a[0],a[1];
h a[0];
s a[0];
s a[0];
h a[1];
s a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
y a[0];
s a[0];
h a[0];
s a[0];
h a[1];
s a[1];
h a[1];
s a[1];
cx a[0],a[1];
h a[0];
s a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
z a[0];
h a[0];
s a[0];
y a[1];
h a[1];
s a[1];
cx a[0],a[1];
h a[0];
s a[0];
h a[1];
s a[1];
cx a[0],a[1];
x a[0];
sdg a[0];
h a[0];
y a[1];
cx a[0],a[1];
z a[0];
s a[0];
z a[1];
cx a[0],a[1];
h a[0];
s a[0];
h a[1];
s a[1];
cx a[0],a[1];
cx a[0],a[1];
y a[0];
s a[0];
h a[0];
s a[0];
y a[1];
s a[1];
h a[1];
s a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
y a[0];
s a[0];
h a[0];
s a[0];
h a[1];
s a[1];
sdg a[1];
h a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
h a[0];
s a[0];
sdg a[0];
h a[0];
h a[0];
y a[1];
s a[1];
z a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
sdg a[0];
h a[0];
sdg a[1];
h a[1];
h a[1];
s a[1];
cx a[0],a[1];
h a[0];
s a[0];
h a[1];
s a[1];
cx a[0],a[1];
y a[0];
h a[0];
s a[0];
x a[1];
cx a[0],a[1];
h a[0];
s a[0];
h a[1];
s a[1];
cx a[0],a[1];
x a[0];
s a[0];
h a[0];
s a[0];
x a[1];
h a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
h a[0];
s a[0];
h a[1];
s a[1];
cx a[0],a[1];
x a[0];
sdg a[0];
h a[0];
y a[1];
s a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
h a[0];
y a[1];
cx a[0],a[1];
z a[0];
h a[0];
z a[1];
h a[1];
s a[1];
cx a[0],a[1];
h a[0];
s a[0];
z a[0];
h a[0];
s a[0];
h a[1];
s a[1];
x a[1];
sdg a[1];
h a[1];
cx a[0],a[1];
h a[0];
s a[0];
z a[0];
sdg a[0];
h a[0];
h a[1];
s a[1];
z a[1];
h a[1];
cx a[0],a[1];
h a[1];
s a[1];
x a[1];
cx a[0],a[1];
x a[0];
s a[0];
sdg a[1];
h a[1];
x a[1];
h a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
h a[1];
s a[1];
cx a[0],a[1];
z a[0];
h a[0];
s a[0];
z a[0];
h a[0];
z a[1];
s a[1];
h a[1];
s a[1];
y a[1];
h a[1];
s a[1];
cx a[0],a[1];
h a[0];
s a[0];
sdg a[1];
h a[1];
z a[1];
sdg a[1];
h a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
h a[1];
s a[1];
cx a[0],a[1];
y a[0];
y a[1];
h a[1];
cx a[0],a[1];
h a[0];
s a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
x a[0];
h a[0];
h a[1];
s a[1];
x a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
y a[0];
s a[0];
h a[0];
s a[0];
h a[1];
s a[1];
cx a[0],a[1];
y a[0];
h a[0];
s a[0];
sdg a[1];
h a[1];
x a[1];
s a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
h a[1];
s a[1];
cx a[0],a[1];
z a[0];
sdg a[0];
h a[0];
z a[1];
s a[1];
h a[1];
s a[1];
cx a[0],a[1];
x a[0];
y a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
h a[1];
s a[1];
cx a[0],a[1];
h a[0];
s a[0];
z a[0];
y a[1];
h a[1];
s a[1];
cx a[0],a[1];
h a[0];
s a[0];
x a[0];
h a[1];
s a[1];
sdg a[1];
h a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
h a[1];
s a[1];
cx a[0],a[1];
h a[0];
x a[0];
z a[1];
h a[1];
s a[1];
x a[1];
cx a[0],a[1];
h a[0];
s a[0];
h a[1];
s a[1];
cx a[0],a[1];
z a[0];
sdg a[0];
h a[0];
z a[1];
cx a[0],a[1];
y a[0];
s a[0];
x a[1];
cx a[0],a[1];
h a[0];
s a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
x a[0];
sdg a[0];
h a[0];
h a[1];
s a[1];
y a[1];
h a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
h a[1];
s a[1];
cx a[0],a[1];
x a[0];
sdg a[0];
h a[0];
x a[1];
s a[1];
h a[1];
s a[1];
cx a[0],a[1];
h a[0];
s a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
h a[1];
s a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
h a[0];
s a[0];
h a[1];
s a[1];
cx a[0],a[1];
x a[0];
h a[0];
s a[0];
s a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
y a[0];
s a[0];
z a[1];
s a[1];
h a[1];
s a[1];
cx a[0],a[1];
y a[0];
sdg a[0];
h a[0];
y a[0];
s a[0];
z a[1];
s a[1];
x a[1];
sdg a[1];
h a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
h a[1];
s a[1];
cx a[0],a[1];
h a[0];
s a[0];
z a[0];
s a[0];
h a[0];
s a[0];
y a[1];
s a[1];
h a[1];
s a[1];
cx a[0],a[1];
h a[0];
s a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
x a[0];
s a[0];
h a[1];
s a[1];
y a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
h a[0];
s a[0];
h a[1];
s a[1];
cx a[0],a[1];
z a[0];
x a[1];
sdg a[1];
h a[1];
cx a[0],a[1];
h a[0];
s a[0];
h a[1];
s a[1];
cx a[0],a[1];
z a[0];
s a[0];
x a[1];
h a[1];
s a[1];
cx a[0],a[1];
h a[0];
s a[0];
y a[0];
sdg a[0];
h a[0];
sdg a[1];
h a[1];
z a[1];
h a[1];
cx a[0],a[1];
h a[0];
s a[0];
h a[1];
s a[1];
cx a[0],a[1];
x a[0];
s a[0];
h a[0];
s a[0];
s a[1];
h a[1];
s a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
h a[0];
h a[1];
s a[1];
y a[1];
s a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
y a[0];
s a[0];
h a[0];
s a[0];
x a[1];
h a[1];
s a[1];
cx a[0],a[1];
h a[0];
s a[0];
h a[0];
sdg a[1];
h a[1];
z a[1];
s a[1];
cx a[0],a[1];
h a[0];
s a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
x a[0];
s a[0];
h a[0];
s a[0];
x a[1];
h a[1];
s a[1];
y a[1];
s a[1];
h a[1];
s a[1];
cx a[0],a[1];
h a[0];
s a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
z a[0];
s a[0];
x a[1];
cx a[0],a[1];
y a[0];
h a[0];
s a[1];
cx a[0],a[1];
z a[0];
z a[1];
h a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
h a[1];
s a[1];
cx a[0],a[1];
z a[0];
h a[0];
s a[0];
h a[0];
y a[1];
s a[1];
x a[1];
h a[1];
s a[1];
cx a[0],a[1];
z a[0];
sdg a[0];
h a[0];
y a[1];
h a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
h a[0];
s a[0];
h a[0];
s a[0];
h a[1];
s a[1];
z a[1];
h a[1];
s a[1];
cx a[0],a[1];
h a[0];
s a[0];
h a[1];
s a[1];
cx a[0],a[1];
y a[0];
x a[1];
h a[1];
s a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
y a[0];
h a[0];
s a[0];
y a[1];
h a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
y a[0];
h a[0];
y a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
x a[0];
h a[0];
x a[1];
h a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
h a[1];
s a[1];
cx a[0],a[1];
z a[0];
h a[0];
h a[1];
s a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
z a[0];
h a[0];
s a[0];
y a[1];
h a[1];
cx a[0],a[1];
h a[0];
s a[0];
z a[0];
s a[0];
h a[0];
s a[0];
x a[1];
sdg a[1];
h a[1];
cx a[0],a[1];
h a[0];
s a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
y a[0];
h a[0];
z a[1];
s a[1];
h a[1];
s a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
z a[0];
h a[0];
sdg a[1];
h a[1];
z a[1];
s a[1];
h a[1];
s a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
y a[0];
h a[0];
s a[0];
x a[1];
s a[1];
h a[1];
s a[1];
cx a[0],a[1];
h a[0];
s a[0];
y a[0];
h a[0];
s a[0];
h a[1];
s a[1];
y a[1];
sdg a[1];
h a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
x a[0];
sdg a[0];
h a[0];
h a[1];
s a[1];
y a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
s a[0];
x a[1];
sdg a[1];
h a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
h a[1];
s a[1];
cx a[0],a[1];
h a[0];
s a[0];
s a[0];
x a[1];
cx a[0],a[1];
s a[0];
h a[0];
s a[0];
h a[1];
s a[1];
sdg a[1];
h a[1];
cx a[0],a[1];
h a[0];
s a[0];
s a[0];
x a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
y a[0];
s a[0];
h a[0];
s a[0];
h a[1];
s a[1];
x a[1];
sdg a[1];
h a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
h a[0];
s a[0];
z a[0];
sdg a[0];
h a[0];
y a[1];
s a[1];
h a[1];
s a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
h a[0];
s a[0];
h a[1];
s a[1];
cx a[0],a[1];
z a[0];
z a[1];
s a[1];
cx a[0],a[1];
h a[0];
s a[0];
s a[0];
h a[1];
s a[1];
y a[1];
h a[1];
cx a[0],a[1];
h a[0];
s a[0];
z a[0];
h a[0];
s a[0];
h a[1];
s a[1];
h a[1];
cx a[0],a[1];
h a[0];
s a[0];
h a[0];
sdg a[1];
h a[1];
y a[1];
s a[1];
h a[1];
s a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
y a[0];
h a[0];
y a[1];
h a[1];
cx a[0],a[1];
x a[0];
h a[0];
s a[0];
h a[1];
s a[1];
y a[1];
h a[1];
s a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
h a[1];
s a[1];
cx a[0],a[1];
y a[0];
y a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
h a[0];
s a[0];
x a[0];
h a[0];
s a[0];
h a[1];
s a[1];
x a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
h a[1];
s a[1];
cx a[0],a[1];
x a[0];
h a[1];
cx a[0],a[1];
h a[0];
s a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
s a[0];
x a[1];
h a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
h a[1];
s a[1];
cx a[0],a[1];
h a[0];
s a[0];
x a[1];
h a[1];
s a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
z a[0];
h a[0];
s a[0];
sdg a[1];
h a[1];
sdg a[1];
h a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
y a[0];
s a[0];
h a[0];
s a[0];
h a[1];
s a[1];
h a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
h a[0];
s a[0];
x a[0];
s a[0];
h a[1];
s a[1];
z a[1];
s a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
h a[1];
s a[1];
y a[1];
cx a[0],a[1];
z a[0];
s a[0];
sdg a[1];
h a[1];
y a[1];
sdg a[1];
h a[1];
cx a[0],a[1];
h a[0];
s a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
z a[0];
sdg a[0];
h a[0];
z a[1];
sdg a[1];
h a[1];
cx a[0],a[1];
h a[0];
s a[0];
h a[1];
s a[1];
cx a[0],a[1];
z a[0];
h a[0];
s a[0];
x a[1];
sdg a[1];
h a[1];
cx a[0],a[1];
h a[0];
s a[0];
h a[1];
s a[1];
cx a[0],a[1];
x a[0];
z a[1];
h a[1];
cx a[0],a[1];
z a[0];
h a[1];
s a[1];
cx a[0],a[1];
h a[0];
s a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
z a[0];
h a[1];
s a[1];
y a[1];
s a[1];
h a[1];
s a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
h a[0];
s a[0];
y a[0];
sdg a[0];
h a[0];
y a[0];
s a[0];
h a[1];
s a[1];
x a[1];
sdg a[1];
h a[1];
cx a[0],a[1];
z a[0];
s a[0];
h a[1];
s a[1];
y a[1];
h a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
h a[1];
s a[1];
cx a[0],a[1];
h a[0];
s a[0];
s a[0];
h a[0];
s a[0];
x a[1];
s a[1];
h a[1];
s a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
z a[0];
sdg a[0];
h a[0];
sdg a[1];
h a[1];
y a[1];
s a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
x a[0];
s a[0];
h a[1];
s a[1];
z a[1];
h a[1];
cx a[0],a[1];
x a[0];
sdg a[0];
h a[0];
sdg a[1];
h a[1];
y a[1];
cx a[0],a[1];
y a[0];
sdg a[0];
h a[0];
h a[0];
h a[1];
s a[1];
x a[1];
sdg a[1];
h a[1];
y a[1];
h a[1];
s a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
z a[0];
sdg a[0];
h a[0];
h a[1];
s a[1];
y a[1];
s a[1];
cx a[0],a[1];
h a[0];
s a[0];
z a[0];
h a[0];
s a[0];
y a[1];
s a[1];
h a[1];
s a[1];
cx a[0],a[1];
h a[0];
s a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
x a[0];
h a[0];
s a[0];
x a[1];
sdg a[1];
h a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
y a[0];
s a[0];
h a[1];
s a[1];
x a[1];
s a[1];
h a[1];
s a[1];
cx a[0],a[1];
h a[0];
s a[0];
z a[0];
h a[0];
z a[1];
h a[1];
s a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
h a[0];
s a[0];
h a[1];
s a[1];
cx a[0],a[1];
s a[0];
x a[1];
h a[1];
s a[1];
cx a[0],a[1];
h a[0];
s a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
z a[0];
s a[0];
h a[0];
s a[0];
h a[1];
s a[1];
y a[1];
h a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
x a[0];
s a[0];
h a[1];
s a[1];
x a[1];
h a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
h a[1];
s a[1];
cx a[0],a[1];
h a[0];
s a[0];
y a[0];
h a[0];
y a[1];
s a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
sdg a[0];
h a[0];
sdg a[1];
h a[1];
x a[1];
h a[1];
s a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
h a[0];
s a[0];
x a[0];
sdg a[0];
h a[0];
h a[1];
s a[1];
x a[1];
h a[1];
cx a[0],a[1];
h a[0];
s a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
x a[0];
s a[0];
x a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
y a[0];
h a[0];
h a[1];
s a[1];
z a[1];
s a[1];
h a[1];
s a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
h a[0];
s a[0];
z a[0];
sdg a[0];
h a[0];
h a[1];
s a[1];
x a[1];
h a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
h a[0];
s a[0];
z a[0];
s a[0];
h a[0];
s a[0];
x a[1];
sdg a[1];
h a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
h a[0];
h a[1];
s a[1];
s a[1];
cx a[0],a[1];
h a[0];
s a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
y a[0];
s a[0];
y a[1];
h a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
z a[0];
h a[0];
s a[0];
sdg a[1];
h a[1];
y a[1];
sdg a[1];
h a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
h a[0];
s a[0];
x a[0];
h a[0];
cx a[0],a[1];
sdg a[0];
h a[0];
h a[1];
s a[1];
cx a[0],a[1];
x a[0];
x a[1];
h a[1];
cx a[0],a[1];
h a[0];
s a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
s a[0];
y a[1];
h a[1];
s a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
h a[0];
s a[0];
h a[0];
s a[0];
x a[1];
sdg a[1];
h a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
h a[0];
s a[0];
x a[0];
h a[0];
s a[0];
h a[1];
s a[1];
z a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
h a[0];
s a[0];
y a[0];
h a[1];
s a[1];
z a[1];
h a[1];
cx a[0],a[1];
z a[0];
x a[1];
h a[1];
s a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
h a[1];
s a[1];
cx a[0],a[1];
h a[0];
s a[0];
z a[0];
sdg a[0];
h a[0];
x a[1];
h a[1];
cx a[0],a[1];
h a[0];
s a[0];
z a[0];
sdg a[1];
h a[1];
z a[1];
s a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
h a[0];
s a[0];
h a[1];
s a[1];
x a[1];
h a[1];
s a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
x a[0];
h a[0];
h a[1];
s a[1];
y a[1];
h a[1];
s a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
z a[0];
s a[0];
h a[0];
s a[0];
h a[1];
s a[1];
h a[1];
cx a[0],a[1];
h a[0];
s a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
y a[0];
s a[0];
h a[0];
s a[0];
h a[1];
s a[1];
x a[1];
sdg a[1];
h a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
x a[0];
sdg a[0];
h a[0];
h a[1];
s a[1];
z a[1];
s a[1];
h a[1];
s a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
h a[1];
s a[1];
cx a[0],a[1];
h a[0];
s a[0];
z a[0];
s a[0];
x a[1];
s a[1];
cx a[0],a[1];
h a[0];
s a[0];
x a[0];
h a[0];
s a[0];
h a[1];
s a[1];
s a[1];
cx a[0],a[1];
z a[0];
s a[0];
h a[1];
s a[1];
x a[1];
s a[1];
h a[1];
s a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
h a[0];
s a[0];
x a[0];
s a[0];
h a[0];
s a[0];
h a[1];
s a[1];
z a[1];
cx a[0],a[1];
z a[0];
s a[0];
sdg a[1];
h a[1];
y a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
y a[0];
h a[1];
s a[1];
sdg a[1];
h a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
z a[0];
s a[0];
x a[1];
s a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
h a[0];
s a[0];
y a[0];
h a[0];
s a[0];
h a[1];
s a[1];
z a[1];
s a[1];
cx a[0],a[1];
y a[0];
s a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
h a[0];
s a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
z a[0];
h a[0];
s a[0];
y a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
h a[1];
s a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
y a[1];
h a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
h a[0];
s a[0];
z a[0];
h a[0];
h a[1];
s a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
z a[0];
s a[0];
x a[1];
sdg a[1];
h a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
h a[1];
s a[1];
cx a[0],a[1];
h a[0];
s a[0];
x a[0];
s a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
h a[1];
s a[1];
cx a[0],a[1];
h a[0];
s a[0];
x a[0];
s a[0];
h a[0];
s a[0];
x a[1];
h a[1];
s a[1];
cx a[0],a[1];
h a[0];
s a[0];
y a[0];
h a[0];
s a[0];
h a[1];
s a[1];
x a[1];
s a[1];
h a[1];
s a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
z a[0];
h a[0];
s a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
s a[0];
h a[0];
s a[0];
x a[1];
s a[1];
h a[1];
s a[1];
cx a[0],a[1];
x a[0];
h a[0];
sdg a[1];
h a[1];
sdg a[1];
h a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
y a[0];
s a[0];
s a[1];
h a[1];
s a[1];
cx a[0],a[1];
z a[0];
sdg a[0];
h a[0];
y a[1];
h a[1];
s a[1];
cx a[0],a[1];
h a[0];
s a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
x a[0];
sdg a[0];
h a[0];
y a[1];
cx a[0],a[1];
h a[0];
s a[0];
x a[0];
h a[0];
s a[0];
h a[1];
s a[1];
z a[1];
s a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
h a[1];
s a[1];
cx a[0],a[1];
h a[0];
s a[0];
x a[0];
h a[0];
h a[1];
cx a[0],a[1];
h a[0];
s a[0];
z a[0];
h a[0];
x a[1];
h a[1];
cx a[0],a[1];
h a[0];
s a[0];
x a[0];
sdg a[0];
h a[0];
h a[1];
s a[1];
y a[1];
sdg a[1];
h a[1];
cx a[0],a[1];
h a[0];
s a[0];
z a[0];
h a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
z a[0];
s a[0];
h a[1];
s a[1];
x a[1];
h a[1];
cx a[0],a[1];
h a[0];
s a[0];
h a[1];
s a[1];
cx a[0],a[1];
x a[0];
sdg a[0];
h a[0];
cx a[0],a[1];
sdg a[0];
h a[0];
z a[0];
s a[0];
h a[0];
s a[0];
sdg a[1];
h a[1];
y a[1];
h a[1];
s a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
x a[1];
s a[1];
cx a[0],a[1];
h a[0];
s a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
y a[0];
h a[0];
s a[1];
h a[1];
s a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
h a[0];
s a[0];
h a[1];
s a[1];
cx a[0],a[1];
h a[0];
x a[1];
h a[1];
s a[1];
cx a[0],a[1];
h a[0];
s a[0];
z a[0];
h a[0];
s a[0];
sdg a[1];
h a[1];
z a[1];
s a[1];
cx a[0],a[1];
h a[0];
s a[0];
x a[0];
sdg a[0];
h a[0];
h a[1];
s a[1];
x a[1];
s a[1];
cx a[0],a[1];
y a[0];
h a[0];
s a[0];
sdg a[1];
h a[1];
h a[1];
y a[1];
s a[1];
h a[1];
s a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
y a[0];
s a[0];
h a[0];
s a[0];
h a[1];
s a[1];
x a[1];
h a[1];
cx a[0],a[1];
x a[0];
sdg a[1];
h a[1];
y a[1];
s a[1];
cx a[0],a[1];
s a[0];
sdg a[1];
h a[1];
z a[1];
s a[1];
cx a[0],a[1];
h a[0];
s a[0];
x a[0];
h a[0];
s a[0];
h a[1];
s a[1];
cx a[0],a[1];
h a[0];
s a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
y a[0];
z a[1];
sdg a[1];
h a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
h a[0];
s a[0];
h a[1];
s a[1];
cx a[0],a[1];
x a[0];
h a[0];
x a[0];
h a[0];
y a[1];
h a[1];
s a[1];
y a[1];
cx a[0],a[1];
h a[0];
s a[0];
h a[1];
s a[1];
cx a[0],a[1];
y a[0];
s a[0];
h a[0];
s a[0];
z a[1];
sdg a[1];
h a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
s a[1];
cx a[0],a[1];
y a[0];
h a[0];
s a[0];
sdg a[1];
h a[1];
x a[1];
cx a[0],a[1];
h a[0];
s a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
x a[0];
s a[0];
h a[0];
s a[0];
x a[0];
s a[0];
h a[0];
s a[0];
h a[1];
y a[1];
sdg a[1];
h a[1];
cx a[0],a[1];
x a[0];
s a[0];
h a[0];
s a[0];
sdg a[1];
h a[1];
s a[1];
cx a[0],a[1];
h a[0];
s a[0];
h a[0];
s a[0];
y a[1];
sdg a[1];
h a[1];
cx a[0],a[1];
h a[0];
s a[0];
x a[0];
h a[0];
s a[0];
h a[1];
s a[1];
z a[1];
h a[1];
s a[1];
cx a[0],a[1];
h a[0];
s a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
z a[0];
h a[0];
s a[0];
x a[1];
s a[1];
h a[1];
s a[1];
cx a[0],a[1];
y a[0];
sdg a[1];
h a[1];
s a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
z a[0];
h a[0];
sdg a[1];
h a[1];
x a[1];
h a[1];
cx a[0],a[1];
h a[0];
s a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
z a[0];
sdg a[0];
h a[0];
z a[1];
s a[1];
h a[1];
s a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
y a[0];
h a[0];
sdg a[1];
h a[1];
s a[1];
h a[1];
s a[1];
cx a[0],a[1];
h a[0];
s a[0];
sdg a[1];
h a[1];
y a[1];
s a[1];
h a[1];
s a[1];
cx a[0],a[1];
h a[0];
s a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
x a[0];
s a[0];
cx a[0],a[1];
h a[0];
s a[0];
x a[0];
s a[0];
h a[0];
s a[0];
sdg a[1];
h a[1];
z a[1];
s a[1];
h a[1];
s a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
h a[0];
s a[0];
x a[0];
s a[0];
x a[0];
h a[0];
h a[1];
s a[1];
h a[1];
s a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
z a[0];
s a[0];
h a[1];
s a[1];
y a[1];
s a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
x a[0];
s a[0];
cx a[0],a[1];
x a[0];
h a[0];
s a[0];
z a[1];
s a[1];
h a[1];
s a[1];
cx a[0],a[1];
h a[0];
s a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
s a[0];
x a[1];
h a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
sdg a[0];
h a[0];
sdg a[1];
h a[1];
z a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
z a[0];
sdg a[0];
h a[0];
z a[1];
h a[1];
cx a[0],a[1];
z a[0];
h a[0];
s a[0];
sdg a[1];
h a[1];
s a[1];
h a[1];
s a[1];
cx a[0],a[1];
x a[0];
h a[0];
sdg a[1];
h a[1];
x a[1];
h a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
h a[1];
s a[1];
cx a[0],a[1];
z a[0];
sdg a[0];
h a[0];
y a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
z a[0];
s a[0];
z a[1];
s a[1];
h a[1];
s a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
y a[0];
sdg a[0];
h a[0];
h a[1];
s a[1];
y a[1];
h a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
h a[0];
s a[0];
z a[0];
s a[0];
h a[1];
s a[1];
y a[1];
h a[1];
s a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
h a[0];
s a[0];
x a[0];
h a[0];
y a[1];
cx a[0],a[1];
h a[0];
s a[0];
y a[0];
h a[0];
sdg a[1];
h a[1];
z a[1];
h a[1];
s a[1];
cx a[0],a[1];
y a[0];
sdg a[0];
h a[0];
sdg a[1];
h a[1];
y a[1];
s a[1];
cx a[0],a[1];
h a[0];
s a[0];
sdg a[0];
h a[0];
h a[1];
s a[1];
z a[1];
h a[1];
cx a[0],a[1];
x a[0];
h a[0];
s a[0];
x a[1];
s a[1];
h a[1];
s a[1];
cx a[0],a[1];
x a[0];
h a[1];
s a[1];
z a[1];
s a[1];
h a[1];
s a[1];
cx a[0],a[1];
h a[0];
s a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
x a[0];
h a[0];
s a[0];
cx a[0],a[1];
sdg a[0];
h a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
x a[0];
s a[0];
h a[1];
cx a[0],a[1];
h a[0];
s a[0];
z a[1];
s a[1];
h a[1];
s a[1];
cx a[0],a[1];
h a[0];
s a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
y a[0];
h a[1];
s a[1];
y a[1];
sdg a[1];
h a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
h a[0];
s a[0];
x a[0];
s a[0];
h a[1];
s a[1];
z a[1];
sdg a[1];
h a[1];
cx a[0],a[1];
y a[0];
h a[0];
h a[1];
s a[1];
s a[1];
cx a[0],a[1];
h a[0];
s a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
z a[0];
s a[0];
h a[0];
s a[0];
y a[1];
sdg a[1];
h a[1];
cx a[0],a[1];
h a[0];
s a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
x a[0];
sdg a[0];
h a[0];
x a[1];
h a[1];
s a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
h a[1];
s a[1];
cx a[0],a[1];
h a[0];
s a[0];
h a[0];
s a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
z a[0];
h a[1];
s a[1];
h a[1];
cx a[0],a[1];
h a[0];
s a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
y a[1];
sdg a[1];
h a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
sdg a[0];
h a[0];
x a[0];
x a[1];
h a[1];
s a[1];
z a[1];
s a[1];
cx a[0],a[1];
y a[0];
h a[0];
s a[0];
h a[1];
s a[1];
x a[1];
s a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
h a[1];
s a[1];
cx a[0],a[1];
h a[0];
s a[0];
z a[0];
x a[1];
h a[1];
cx a[0],a[1];
s a[0];
x a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
h a[0];
s a[0];
y a[0];
s a[0];
h a[1];
s a[1];
y a[1];
s a[1];
h a[1];
s a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
x a[0];
s a[0];
h a[0];
s a[0];
h a[1];
s a[1];
x a[1];
h a[1];
cx a[0],a[1];
h a[0];
s a[0];
x a[0];
h a[0];
s a[0];
z a[1];
s a[1];
h a[1];
s a[1];
cx a[0],a[1];
h a[0];
s a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
z a[0];
s a[0];
h a[1];
s a[1];
x a[1];
cx a[0],a[1];
h a[0];
s a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
y a[0];
h a[0];
z a[1];
h a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
y a[0];
h a[1];
s a[1];
sdg a[1];
h a[1];
cx a[0],a[1];
h a[0];
s a[0];
h a[1];
s a[1];
cx a[0],a[1];
h a[0];
s a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
s a[0];
z a[0];
s a[0];
h a[0];
s a[0];
x a[1];
s a[1];
x a[1];
s a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
h a[0];
s a[0];
h a[1];
s a[1];
h a[1];
s a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
y a[0];
h a[0];
h a[1];
s a[1];
y a[1];
s a[1];
h a[1];
s a[1];
cx a[0],a[1];
h a[0];
s a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
z a[0];
s a[0];
h a[0];
s a[0];
x a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
x a[0];
y a[1];
h a[1];
s a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
x a[0];
s a[0];
y a[1];
s a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
y a[0];
h a[1];
s a[1];
s a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
z a[0];
sdg a[0];
h a[0];
x a[1];
h a[1];
s a[1];
cx a[0],a[1];
h a[0];
s a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
z a[0];
h a[0];
s a[0];
x a[1];
h a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
h a[0];
s a[0];
z a[0];
s a[0];
h a[1];
s a[1];
x a[1];
s a[1];
h a[1];
s a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
h a[1];
s a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
z a[1];
h a[1];
cx a[0],a[1];
h a[0];
s a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
x a[0];
s a[0];
h a[1];
s a[1];
z a[1];
s a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
h a[0];
h a[1];
s a[1];
x a[1];
s a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
h a[1];
s a[1];
y a[1];
h a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
z a[0];
sdg a[0];
h a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
h a[0];
s a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
z a[0];
s a[0];
h a[0];
s a[0];
h a[1];
s a[1];
z a[1];
h a[1];
s a[1];
cx a[0],a[1];
h a[0];
s a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
z a[0];
x a[1];
h a[1];
cx a[0],a[1];
h a[0];
s a[0];
x a[0];
h a[0];
sdg a[1];
h a[1];
z a[1];
h a[1];
s a[1];
cx a[0],a[1];
z a[0];
sdg a[0];
h a[0];
sdg a[1];
h a[1];
y a[1];
s a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
z a[0];
sdg a[1];
h a[1];
x a[1];
h a[1];
cx a[0],a[1];
h a[0];
s a[0];
s a[0];
z a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
x a[0];
s a[0];
h a[0];
s a[0];
h a[1];
s a[1];
y a[1];
h a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
h a[1];
s a[1];
cx a[0],a[1];
z a[0];
y a[1];
cx a[0],a[1];
h a[0];
s a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
h a[1];
s a[1];
s a[1];
h a[1];
s a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
x a[0];
h a[1];
s a[1];
s a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
h a[1];
s a[1];
cx a[0],a[1];
s a[0];
x a[1];
sdg a[1];
h a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
x a[0];
s a[0];
sdg a[1];
h a[1];
sdg a[1];
h a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
z a[0];
h a[0];
sdg a[1];
h a[1];
x a[1];
h a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
y a[0];
h a[0];
s a[0];
h a[1];
s a[1];
s a[1];
cx a[0],a[1];
h a[0];
s a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
y a[0];
h a[0];
h a[1];
s a[1];
z a[1];
h a[1];
cx a[0],a[1];
x a[0];
sdg a[1];
h a[1];
x a[1];
h a[1];
s a[1];
cx a[0],a[1];
h a[0];
s a[0];
h a[1];
s a[1];
cx a[0],a[1];
y a[0];
s a[0];
h a[0];
s a[0];
z a[1];
h a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
y a[0];
s a[0];
x a[1];
s a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
h a[0];
s a[0];
h a[1];
s a[1];
cx a[0],a[1];
x a[0];
h a[0];
h a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
z a[0];
h a[1];
s a[1];
y a[1];
h a[1];
s a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
x a[0];
h a[0];
h a[1];
s a[1];
x a[1];
h a[1];
s a[1];
cx a[0],a[1];
h a[0];
s a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
z a[0];
z a[1];
s a[1];
h a[1];
s a[1];
cx a[0],a[1];
y a[0];
h a[0];
s a[0];
y a[1];
s a[1];
h a[1];
s a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
h a[1];
s a[1];
cx a[0],a[1];
z a[0];
h a[0];
s a[0];
y a[1];
h a[1];
s a[1];
cx a[0],a[1];
h a[0];
s a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
y a[0];
s a[0];
z a[1];
sdg a[1];
h a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
h a[0];
s a[0];
h a[1];
s a[1];
cx a[0],a[1];
y a[0];
y a[1];
sdg a[1];
h a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
h a[0];
s a[0];
y a[0];
h a[0];
z a[1];
h a[1];
s a[1];
cx a[0],a[1];
h a[0];
s a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
z a[0];
sdg a[0];
h a[0];
h a[1];
s a[1];
z a[1];
s a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
h a[1];
s a[1];
cx a[0],a[1];
x a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
h a[0];
s a[0];
h a[1];
s a[1];
cx a[0],a[1];
z a[0];
h a[0];
h a[1];
cx a[0],a[1];
z a[0];
h a[0];
cx a[0],a[1];
sdg a[0];
h a[0];
h a[1];
s a[1];
cx a[0],a[1];
h a[0];
s a[0];
y a[0];
h a[0];
s a[0];
x a[1];
h a[1];
cx a[0],a[1];
h a[0];
s a[0];
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
