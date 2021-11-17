OPENQASM 2.0;
include "qelib1.inc";
gate circuit-7023 q0,q1,q2,q3,q4 { s q1; h q2; cx q0,q3; cx q2,q0; h q2; s q3; s q4; cx q4,q1; cx q1,q0; cx q0,q4; y q0; s q0; h q0; s q1; swap q4,q1; s q1; h q1; cx q4,q2; h q2; rz(1.418176) q2; cx q2,q1; cx q4,q3; h q3; rz(1.2241042) q3; rz(3.9601559) q4; cx q3,q4; h q3; cx q1,q3; s q1; cx q3,q2; h q2; cx q2,q1; s q1; h q1; s q1; s q2; h q2; z q2; h q2; swap q2,q0; cx q2,q0; h q0; rz(3.2450776) q0; y q3; cx q3,q2; rz(2.7590799) q3; s q4; cx q4,q1; }
qreg q4[5];
qreg q5[1];
rx(1.979546) q4[0];
ry(1.0647481) q4[0];
rz(1.6637346) q4[0];
rx(1.952188) q4[1];
ry(1.1686433) q4[1];
rz(1.9456426) q4[1];
rx(1.9665311) q4[2];
ry(0.83558878) q4[2];
rz(0.5501028) q4[2];
rx(1.3332606) q4[3];
ry(1.0030967) q4[3];
rz(0.027446276) q4[3];
rx(0.063510311) q4[4];
ry(0.13124358) q4[4];
rz(0.81823185) q4[4];
h q5[0];
barrier q4[0],q4[1],q4[2],q4[3],q4[4],q5[0];
h q4[4];
cx q5[0],q4[4];
h q4[4];
cx q5[0],q4[3];
cx q5[0],q4[2];
h q4[1];
cx q5[0],q4[1];
h q4[1];
cx q5[0],q4[0];
barrier q4[0],q4[1],q4[2],q4[3],q4[4],q5[0];
circuit-7023 q4[0],q4[1],q4[2],q4[3],q4[4];
barrier q4[0],q4[1],q4[2],q4[3],q4[4],q5[0];
h q4[4];
cx q5[0],q4[4];
h q4[4];
h q4[3];
cx q5[0],q4[3];
h q4[3];
sdg q4[2];
cx q5[0],q4[2];
s q4[2];
cx q5[0],q4[1];
h q4[0];
cx q5[0],q4[0];
h q4[0];
h q5[0];
