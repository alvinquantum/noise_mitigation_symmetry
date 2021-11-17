OPENQASM 2.0;
include "qelib1.inc";
gate circuit-7015 q0,q1,q2,q3,q4 { h q0; rz(2.4882911) q1; swap q1,q0; cx q1,q2; rz(3.1800398) q1; rz(5.5971384) q2; h q3; rz(3.3829821) q3; rz(1.2399226) q4; s q4; }
qreg q4[5];
qreg q5[1];
rx(0.10313029) q4[0];
ry(1.0309844) q4[0];
rz(1.8640753) q4[0];
rx(0.89883394) q4[1];
ry(1.363562) q4[1];
rz(0.80342813) q4[1];
rx(1.3313001) q4[2];
ry(0.34549) q4[2];
rz(0.39351976) q4[2];
rx(0.70630525) q4[3];
ry(1.0279345) q4[3];
rz(0.82394193) q4[3];
rx(1.1784369) q4[4];
ry(0.7030592) q4[4];
rz(0.7352737) q4[4];
h q5[0];
barrier q4[0],q4[1],q4[2],q4[3],q4[4],q5[0];
h q4[4];
cx q5[0],q4[4];
h q4[4];
cx q5[0],q4[3];
h q4[2];
cx q5[0],q4[2];
h q4[2];
h q4[1];
cx q5[0],q4[1];
h q4[1];
barrier q4[0],q4[1],q4[2],q4[3],q4[4],q5[0];
circuit-7015 q4[0],q4[1],q4[2],q4[3],q4[4];
barrier q4[0],q4[1],q4[2],q4[3],q4[4],q5[0];
h q4[4];
cx q5[0],q4[4];
h q4[4];
h q4[3];
cx q5[0],q4[3];
h q4[3];
h q4[2];
cx q5[0],q4[2];
h q4[2];
h q4[1];
cx q5[0],q4[1];
h q4[1];
h q4[0];
cx q5[0],q4[0];
h q4[0];
h q5[0];
