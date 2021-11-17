OPENQASM 2.0;
include "qelib1.inc";
gate circuit-7035 q0,q1,q2,q3,q4 { s q0; h q0; s q0; s q1; h q1; h q2; s q2; h q3; s q4; cx q0,q4; cx q1,q0; cx q0,q3; s q1; h q1; s q1; h q3; cx q3,q0; rz(5.8662516) q0; y q0; s q0; h q0; h q3; swap q3,q1; cx q1,q3; h q1; h q4; cx q4,q2; cx q2,q3; h q2; cx q3,q4; x q3; h q3; cx q3,q0; s q4; h q4; s q4; cx q4,q2; cx q2,q1; cx q1,q4; x q1; s q1; s q2; h q2; s q2; h q4; swap q4,q2; cx q4,q2; h q2; z q2; x q4; s q4; h q4; s q4; cx q2,q4; h q2; cx q4,q1; h q1; cx q1,q4; cx q0,q4; s q0; h q1; cx q4,q3; s q3; cx q0,q3; cx q1,q0; h q1; s q1; cx q2,q0; z q0; s q0; s q2; rz(1.405645) q3; cx q3,q1; cx q1,q2; s q1; cx q2,q3; y q2; h q2; rz(5.9379293) q3; h q3; swap q3,q1; cx q1,q3; s q1; h q1; s q1; h q1; swap q1,q0; cx q1,q2; h q2; h q3; cx q1,q3; cx q1,q4; cx q1,q0; x q1; s q4; h q4; swap q4,q0; cx q4,q3; cx q2,q4; cx q0,q4; s q0; h q0; s q2; cx q2,q0; s q3; h q3; s q3; cx q0,q3; h q0; cx q3,q2; swap q2,q0; rz(0.28253254) q2; cx q2,q0; h q0; cx q0,q2; h q0; z q0; swap q1,q0; s q2; cx q2,q0; s q0; cx q1,q2; h q2; rz(5.9080267) q2; y q3; h q3; s q3; z q4; h q4; s q4; cx q4,q3; }
qreg q4[5];
qreg q5[1];
rx(1.1862852) q4[0];
ry(0.27395057) q4[0];
rz(0.77286101) q4[0];
rx(0.59008904) q4[1];
ry(0.51854368) q4[1];
rz(1.6466453) q4[1];
rx(1.1974682) q4[2];
ry(1.8505441) q4[2];
rz(0.6247622) q4[2];
rx(1.1282734) q4[3];
ry(1.0204148) q4[3];
rz(0.067604326) q4[3];
rx(0.37803924) q4[4];
ry(1.0291236) q4[4];
rz(1.4004193) q4[4];
h q5[0];
barrier q4[0],q4[1],q4[2],q4[3],q4[4],q5[0];
cx q5[0],q4[3];
h q4[2];
cx q5[0],q4[2];
h q4[2];
cx q5[0],q4[1];
barrier q4[0],q4[1],q4[2],q4[3],q4[4],q5[0];
circuit-7035 q4[0],q4[1],q4[2],q4[3],q4[4];
barrier q4[0],q4[1],q4[2],q4[3],q4[4],q5[0];
sdg q4[4];
cx q5[0],q4[4];
s q4[4];
cx q5[0],q4[3];
h q4[2];
cx q5[0],q4[2];
h q4[2];
h q4[1];
cx q5[0],q4[1];
h q4[1];
cx q5[0],q4[0];
h q5[0];
