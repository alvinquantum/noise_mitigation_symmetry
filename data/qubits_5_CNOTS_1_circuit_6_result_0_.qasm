OPENQASM 2.0;
include "qelib1.inc";
gate circuit-10550 q0,q1,q2,q3,q4 { s q0; s q2; cx q1,q2; rz(5.5800263) q1; rz(3.96665) q2; rz(4.6531782) q3; s q3; rz(2.9298605) q4; s q4; h q4; rz(2.0556308) q4; }
qreg q6[5];
qreg q7[1];
rx(0.24054862) q6[0];
ry(1.8356248) q6[0];
rz(1.3050662) q6[0];
rx(0.6709765) q6[1];
ry(0.02280495) q6[1];
rz(0.71107132) q6[1];
rx(1.7223365) q6[2];
ry(0.84720111) q6[2];
rz(1.7313364) q6[2];
rx(0.57529256) q6[3];
ry(0.82657969) q6[3];
rz(1.5407623) q6[3];
rx(1.0815658) q6[4];
ry(1.1780175) q6[4];
rz(0.1777339) q6[4];
h q7[0];
barrier q6[0],q6[1],q6[2],q6[3],q6[4],q7[0];
h q6[3];
cx q7[0],q6[3];
h q6[3];
h q6[2];
cx q7[0],q6[2];
h q6[2];
cx q7[0],q6[0];
barrier q6[0],q6[1],q6[2],q6[3],q6[4],q7[0];
circuit-10550 q6[0],q6[1],q6[2],q6[3],q6[4];
barrier q6[0],q6[1],q6[2],q6[3],q6[4],q7[0];
h q6[3];
cx q7[0],q6[3];
h q6[3];
h q6[2];
cx q7[0],q6[2];
h q6[2];
h q6[1];
cx q7[0],q6[1];
h q6[1];
sdg q6[0];
cx q7[0],q6[0];
s q6[0];
h q7[0];
