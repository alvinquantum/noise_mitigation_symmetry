OPENQASM 2.0;
include "qelib1.inc";
qreg a[2];
qreg b[1];
id a[0];
id a[1];
id b[0];
rx(0.67690137) a[0];
ry(0.52222697) a[0];
rz(0.77315361) a[0];
rx(1.8504465) a[1];
ry(1.7822107) a[1];
rz(1.6550967) a[1];
barrier a[0],a[1],b[0];
id b[0];
h b[0];
id a[0];
id a[1];
h a[1];
cx b[0],a[1];
h a[1];
barrier a[0],a[1],b[0];
id b[0];
id a[0];
y a[0];
h a[0];
id a[1];
z a[1];
sdg a[1];
h a[1];
cx a[0],a[1];
h a[0];
s a[0];
h a[1];
s a[1];
x a[1];
s a[1];
cx a[0],a[1];
h a[0];
s a[0];
z a[0];
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
y a[0];
h a[0];
s a[0];
h a[1];
s a[1];
z a[1];
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
x a[1];
cx a[0],a[1];
h a[0];
s a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
x a[0];
h a[0];
z a[1];
s a[1];
cx a[0],a[1];
h a[0];
s a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
x a[0];
y a[1];
sdg a[1];
h a[1];
cx a[0],a[1];
z a[0];
z a[1];
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
x a[1];
s a[1];
cx a[0],a[1];
h a[0];
s a[0];
z a[0];
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
h a[0];
s a[0];
sdg a[0];
h a[0];
z a[1];
s a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
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
