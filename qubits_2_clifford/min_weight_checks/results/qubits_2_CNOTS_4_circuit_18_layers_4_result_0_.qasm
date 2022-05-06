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
rx(0.30760979) a[0];
ry(0.052590709) a[0];
rz(1.2307198) a[0];
rx(1.8702539) a[1];
ry(1.1737499) a[1];
rz(0.40768431) a[1];
barrier a[0],a[1],b[0],b[1],b[2],b[3];
id a[0];
sdg a[0];
id a[1];
u3(pi,5.8991791,2.7575864) a[1];
id b[0];
h b[0];
id b[1];
h b[1];
u3(1.0592554e-13,3.288277,2.8116442) b[1];
id b[2];
h b[2];
u3(1.0592554e-13,3.288277,2.8116442) b[2];
id b[3];
h b[3];
u3(1.0592554e-13,3.288277,2.8116442) b[3];
cx b[3],a[1];
u3(pi,0.38400625,3.5255989) a[1];
h a[1];
u3(pi,5.8991791,2.7575864) a[1];
cx b[1],a[1];
u3(pi,0.38400625,3.5255989) a[1];
h a[1];
cx b[0],a[1];
u3(2*pi,5.4995268,4.1085153) b[1];
u3(2*pi,5.4995268,4.1085153) b[3];
cx b[3],a[0];
s a[0];
u3(pi,5.8991791,2.7575864) a[0];
cx b[2],a[0];
u3(pi,0.38400625,3.5255989) a[0];
cx b[1],a[0];
u3(2*pi,5.4995268,4.1085153) b[2];
barrier a[0],a[1],b[0],b[1],b[2],b[3];
id a[0];
z a[0];
s a[0];
h a[0];
s a[0];
id a[1];
y a[1];
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
h a[0];
x a[1];
sdg a[1];
h a[1];
cx a[0],a[1];
h a[0];
s a[0];
y a[0];
s a[0];
s a[1];
cx a[0],a[1];
h a[0];
s a[0];
h a[1];
s a[1];
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
