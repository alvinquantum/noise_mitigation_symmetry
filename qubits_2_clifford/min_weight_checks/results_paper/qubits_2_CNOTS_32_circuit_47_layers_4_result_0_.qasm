OPENQASM 2.0;
include "qelib1.inc";
qreg a[2];
qreg b[4];
id a[0];
id a[1];
id b[0];
id b[1];
id b[2];
id b[3];
rx(1.4658185) a[0];
ry(0.44047866) a[0];
rz(0.86864256) a[0];
rx(1.3261465) a[1];
ry(0.59136529) a[1];
rz(0.18543938) a[1];
barrier a[0],a[1],b[0],b[1],b[2],b[3];
id a[0];
sdg a[0];
id a[1];
id b[0];
h b[0];
id b[1];
h b[1];
cx b[1],a[0];
s a[0];
h a[0];
cx b[0],a[0];
h a[0];
id b[2];
h b[2];
id b[3];
h b[3];
cx b[3],a[1];
h a[1];
cx b[2],a[1];
h a[1];
barrier a[0],a[1],b[0],b[1],b[2],b[3];
id a[0];
h a[0];
id a[1];
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
h a[0];
s a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
h a[1];
s a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
x a[1];
sdg a[1];
h a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
y a[0];
sdg a[0];
h a[0];
h a[1];
s a[1];
y a[1];
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
x a[1];
sdg a[1];
h a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
h a[1];
s a[1];
cx a[0],a[1];
x a[0];
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
z a[0];
s a[0];
h a[1];
s a[1];
y a[1];
sdg a[1];
h a[1];
cx a[0],a[1];
h a[0];
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
h a[1];
s a[1];
cx a[0],a[1];
h a[0];
s a[0];
x a[0];
h a[0];
z a[1];
h a[1];
cx a[0],a[1];
h a[0];
s a[0];
s a[0];
h a[0];
s a[0];
x a[1];
h a[1];
s a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
sdg a[0];
h a[0];
h a[1];
s a[1];
s a[1];
h a[1];
s a[1];
cx a[0],a[1];
h a[0];
s a[0];
x a[0];
z a[1];
cx a[0],a[1];
sdg a[0];
h a[0];
h a[0];
s a[0];
sdg a[1];
h a[1];
x a[1];
cx a[0],a[1];
x a[0];
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
s a[0];
h a[0];
s a[0];
h a[1];
s a[1];
cx a[0],a[1];
h a[0];
s a[0];
x a[0];
s a[0];
z a[0];
s a[0];
h a[0];
s a[0];
h a[1];
s a[1];
cx a[0],a[1];
h a[0];
s a[0];
y a[0];
h a[1];
s a[1];
x a[1];
h a[1];
s a[1];
cx a[0],a[1];
h a[0];
s a[0];
x a[0];
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
x a[0];
h a[0];
z a[0];
s a[0];
h a[1];
s a[1];
y a[1];
s a[1];
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
id b[0];
id b[1];
id b[2];
id b[3];
barrier a[0],a[1],b[0],b[1],b[2],b[3];
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
h a[0];
h b[2];
id b[3];
cx b[3],a[0];
h a[0];
h b[3];
barrier a[0],a[1],b[0],b[1],b[2],b[3];
id a[0];
id a[1];
id b[0];
id b[1];
id b[2];
id b[3];
