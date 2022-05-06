OPENQASM 2.0;
include "qelib1.inc";
qreg a[2];
qreg b[1];
id a[0];
id a[1];
id b[0];
rx(1.8881762) a[0];
ry(0.24926315) a[0];
rz(1.8286747) a[0];
rx(1.277725) a[1];
ry(1.3999925) a[1];
rz(1.5080761) a[1];
barrier a[0],a[1],b[0];
id b[0];
h b[0];
id a[0];
sdg a[0];
id a[1];
cx b[0],a[1];
cx b[0],a[0];
s a[0];
barrier a[0],a[1],b[0];
id b[0];
id a[0];
z a[0];
h a[0];
id a[1];
z a[1];
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
h a[0];
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
z a[1];
s a[1];
cx a[0],a[1];
h a[0];
s a[0];
h a[1];
s a[1];
z a[1];
h a[1];
cx a[0],a[1];
x a[0];
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
h a[1];
s a[1];
cx a[0],a[1];
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
sdg a[0];
h a[0];
y a[1];
h a[1];
cx a[0],a[1];
h a[0];
s a[0];
z a[0];
sdg a[0];
h a[0];
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
y a[0];
s a[0];
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
x a[1];
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
cx a[0],a[1];
h a[0];
s a[0];
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
sdg a[1];
h a[1];
cx a[0],a[1];
z a[0];
h a[0];
z a[0];
h a[1];
s a[1];
y a[1];
h a[1];
y a[1];
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
sdg a[0];
h a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
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
y a[0];
sdg a[0];
h a[0];
z a[1];
sdg a[1];
h a[1];
z a[1];
sdg a[1];
h a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
x a[0];
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
h a[1];
s a[1];
cx a[0],a[1];
x a[0];
s a[0];
y a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
h a[1];
s a[1];
cx a[0],a[1];
z a[0];
x a[1];
h a[1];
s a[1];
cx a[0],a[1];
h a[0];
s a[0];
z a[0];
s a[0];
h a[0];
s a[0];
h a[1];
s a[1];
x a[1];
cx a[0],a[1];
h a[0];
s a[0];
y a[0];
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
z a[0];
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
s a[0];
h a[0];
s a[0];
x a[1];
s a[1];
cx a[0],a[1];
h a[0];
s a[0];
y a[0];
s a[0];
h a[0];
s a[0];
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
z a[0];
z a[1];
h a[1];
s a[1];
cx a[0],a[1];
y a[0];
sdg a[0];
h a[0];
x a[0];
h a[0];
s a[0];
z a[1];
s a[1];
h a[1];
s a[1];
x a[1];
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
s a[0];
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
h a[1];
s a[1];
cx a[0],a[1];
z a[0];
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
x a[0];
h a[0];
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
z a[0];
cx a[0],a[1];
sdg a[0];
h a[0];
h a[1];
s a[1];
cx a[0],a[1];
z a[0];
s a[0];
h a[0];
s a[0];
z a[1];
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
x a[1];
h a[1];
s a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
h a[1];
s a[1];
cx a[0],a[1];
z a[0];
x a[1];
h a[1];
s a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
h a[1];
s a[1];
z a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
x a[0];
s a[0];
sdg a[1];
h a[1];
h a[1];
cx a[0],a[1];
h a[0];
s a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
z a[0];
h a[0];
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
y a[0];
h a[0];
h a[1];
s a[1];
x a[1];
sdg a[1];
h a[1];
cx a[0],a[1];
h a[0];
s a[0];
z a[0];
h a[0];
sdg a[1];
h a[1];
y a[1];
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
z a[1];
s a[1];
cx a[0],a[1];
h a[0];
s a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
y a[1];
s a[1];
h a[1];
s a[1];
cx a[0],a[1];
y a[0];
s a[0];
h a[0];
s a[0];
y a[1];
sdg a[1];
h a[1];
cx a[0],a[1];
h a[0];
s a[0];
s a[0];
x a[0];
z a[1];
s a[1];
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
z a[0];
h a[0];
x a[1];
h a[1];
cx a[0],a[1];
h a[0];
s a[0];
z a[0];
s a[0];
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
z a[0];
s a[0];
y a[1];
cx a[0],a[1];
y a[0];
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
z a[0];
sdg a[0];
h a[0];
s a[1];
h a[1];
s a[1];
cx a[0],a[1];
h a[0];
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
h a[0];
s a[0];
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
h a[1];
s a[1];
cx a[0],a[1];
s a[0];
h a[0];
s a[0];
y a[1];
s a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
x a[0];
sdg a[0];
h a[0];
h a[1];
s a[1];
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
z a[0];
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
s a[0];
sdg a[1];
h a[1];
y a[1];
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
h a[0];
s a[0];
y a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
x a[0];
s a[0];
h a[0];
s a[0];
x a[1];
h a[1];
cx a[0],a[1];
x a[0];
s a[0];
s a[1];
h a[1];
s a[1];
cx a[0],a[1];
h a[0];
s a[0];
h a[1];
s a[1];
cx a[0],a[1];
z a[0];
h a[0];
y a[1];
sdg a[1];
h a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
h a[0];
sdg a[1];
h a[1];
x a[1];
s a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
z a[0];
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
sdg a[1];
h a[1];
cx a[0],a[1];
h a[0];
s a[0];
s a[0];
h a[1];
s a[1];
x a[1];
sdg a[1];
h a[1];
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
s a[0];
h a[0];
s a[0];
x a[1];
h a[1];
cx a[0],a[1];
x a[0];
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
z a[1];
h a[1];
x a[1];
h a[1];
cx a[0],a[1];
y a[0];
s a[0];
h a[0];
s a[0];
z a[1];
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
h a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
sdg a[1];
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
h a[0];
s a[0];
h a[1];
s a[1];
cx a[0],a[1];
s a[0];
x a[1];
h a[1];
cx a[0],a[1];
z a[0];
sdg a[0];
h a[0];
z a[1];
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
h a[0];
s a[0];
z a[0];
s a[0];
h a[0];
s a[0];
z a[1];
s a[1];
cx a[0],a[1];
z a[0];
h a[1];
s a[1];
y a[1];
h a[1];
cx a[0],a[1];
h a[0];
s a[0];
sdg a[1];
h a[1];
cx a[0],a[1];
y a[0];
sdg a[0];
h a[0];
h a[1];
s a[1];
x a[1];
h a[1];
cx a[0],a[1];
h a[0];
s a[0];
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
