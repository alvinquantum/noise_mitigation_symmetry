OPENQASM 2.0;
include "qelib1.inc";
gate circuit-9 q0,q1,q2,q3,q4 { h q0; s q1; h q1; s q2; h q2; cx q1,q2; s q1; h q3; cx q3,q0; cx q0,q2; s q0; h q0; s q0; cx q2,q3; y q2; s q2; swap q3,q0; cx q0,q3; h q0; rz(1.0923141) q0; cx q3,q1; h q1; cx q1,q3; swap q1,q0; rz(5.7547059) q1; cx q0,q1; h q0; y q3; h q3; cx q1,q4; h q4; cx q4,q1; z q1; s q1; h q1; h q4; swap q4,q0; cx q0,q4; s q0; h q0; y q0; h q0; s q0; cx q0,q2; cx q1,q0; s q1; h q1; s q1; swap q3,q1; cx q2,q3; h q2; rz(0.94165346) q4; z q4; h q4; cx q4,q0; s q0; h q4; cx q4,q3; cx q3,q1; h q1; cx q1,q3; s q1; z q3; s q3; h q3; cx q3,q0; h q4; cx q4,q2; cx q2,q1; cx q1,q4; z q1; h q1; h q2; rz(0.95575014) q2; h q4; cx q4,q2; h q4; s q4; x q4; cx q4,q2; cx q2,q1; h q1; cx q1,q2; cx q0,q2; cx q2,q3; y q2; h q2; s q3; cx q3,q0; cx q0,q1; s q0; cx q0,q4; cx q1,q3; z q1; s q1; h q1; s q3; h q3; h q4; cx q4,q0; s q0; h q0; h q4; cx q3,q4; z q3; h q3; cx q0,q3; rz(1.8878859) q3; s q4; h q4; z q4; h q4; }
qreg q0[5];
qreg q1[1];
rx(1.3879512) q0[0];
ry(0.12881282) q0[0];
rz(0.98074668) q0[0];
rx(1.5695849) q0[1];
ry(1.0859509) q0[1];
rz(0.79823084) q0[1];
rx(1.3042811) q0[2];
ry(1.299944) q0[2];
rz(1.9087051) q0[2];
rx(1.0923028) q0[3];
ry(1.4634389) q0[3];
rz(1.3235965) q0[3];
rx(1.171262) q0[4];
ry(1.1941733) q0[4];
rz(1.0417292) q0[4];
h q1[0];
barrier q0[0],q0[1],q0[2],q0[3],q0[4],q1[0];
sdg q0[3];
cx q1[0],q0[3];
s q0[3];
cx q1[0],q0[2];
h q0[1];
cx q1[0],q0[1];
h q0[1];
h q0[0];
cx q1[0],q0[0];
h q0[0];
barrier q0[0],q0[1],q0[2],q0[3],q0[4],q1[0];
circuit-9 q0[0],q0[1],q0[2],q0[3],q0[4];
barrier q0[0],q0[1],q0[2],q0[3],q0[4],q1[0];
cx q1[0],q0[4];
h q0[3];
cx q1[0],q0[3];
h q0[3];
h q0[2];
cx q1[0],q0[2];
h q0[2];
sdg q0[1];
cx q1[0],q0[1];
s q0[1];
sdg q0[0];
cx q1[0],q0[0];
s q0[0];
h q1[0];
