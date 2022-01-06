OPENQASM 2.0;
include "qelib1.inc";
qreg q186[5];
qreg q187[1];
h q186[4];
cx q187[0],q186[4];
h q186[4];
cx q187[0],q186[3];
cx q187[0],q186[1];
sdg q186[0];
cx q187[0],q186[0];
s q186[0];
barrier q186[0],q186[1],q186[2],q186[3],q186[4],q187[0];
s q186[0];
s q186[2];
h q186[3];
rz(2.9693335) q186[4];
s q186[4];
h q186[0];
rz(5.3177183) q186[2];
cx q186[2],q186[4];
s q186[0];
rz(0.58840822) q186[4];
h q186[4];
rz(3.1988129) q186[2];
rz(2.6568052) q186[0];
barrier q186[0],q186[1],q186[2],q186[3],q186[4],q187[0];
u3(1.0592554e-13,3.288277,2.8116442) q187[0];
u3(pi,5.8991791,2.7575864) q186[4];
cx q187[0],q186[4];
u3(2*pi,5.4995268,4.1085153) q187[0];
u3(pi,0.38400625,3.5255989) q186[4];
h q186[3];
cx q187[0],q186[3];
h q186[3];
h q186[2];
cx q187[0],q186[2];
h q186[2];
cx q187[0],q186[1];
h q186[0];
cx q187[0],q186[0];
h q186[0];
