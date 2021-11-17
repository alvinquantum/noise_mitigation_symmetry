OPENQASM 2.0;
include "qelib1.inc";
gate circuit-9 q0,q1,q2,q3,q4 { s q0; h q0; swap q1,q0; h q3; s q3; cx q3,q2; h q4; rz(3.244989) q4; cx q4,q1; cx q1,q0; h q0; rz(2.9491007) q0; cx q0,q1; s q0; h q0; s q0; rz(3.1792464) q0; cx q2,q1; cx q1,q3; s q2; h q2; h q3; cx q3,q0; s q3; cx q4,q2; cx q2,q0; cx q0,q4; s q0; h q0; s q0; s q2; h q2; s q4; h q4; cx q4,q2; rz(0.55161358) q2; cx q2,q3; h q3; cx q3,q2; y q2; s q3; h q3; s q4; swap q4,q3; cx q4,q3; s q3; h q3; y q3; x q4; s q4; swap q4,q0; cx q4,q1; h q1; cx q1,q2; h q2; cx q4,q0; s q0; rz(5.8365237) q0; h q0; cx q4,q3; h q3; cx q3,q4; s q3; h q3; cx q3,q0; cx q0,q1; h q1; cx q1,q0; y q0; h q1; cx q2,q1; h q3; s q3; cx q1,q3; h q1; cx q3,q2; h q2; swap q2,q1; y q3; }
qreg q0[5];
qreg q1[1];
rx(0.58080466) q0[0];
ry(1.2761056) q0[0];
rz(0.66779021) q0[0];
rx(0.64120898) q0[1];
ry(0.12872149) q0[1];
rz(0.88644464) q0[1];
rx(1.149386) q0[2];
ry(1.5768936) q0[2];
rz(0.88821086) q0[2];
rx(0.079785089) q0[3];
ry(0.44684944) q0[3];
rz(0.44865114) q0[3];
rx(0.34900348) q0[4];
ry(1.3339915) q0[4];
rz(0.47157886) q0[4];
h q1[0];
barrier q0[0],q0[1],q0[2],q0[3],q0[4],q1[0];
sdg q0[3];
cx q1[0],q0[3];
s q0[3];
h q0[2];
cx q1[0],q0[2];
h q0[2];
cx q1[0],q0[0];
barrier q0[0],q0[1],q0[2],q0[3],q0[4],q1[0];
circuit-9 q0[0],q0[1],q0[2],q0[3],q0[4];
barrier q0[0],q0[1],q0[2],q0[3],q0[4],q1[0];
h q0[4];
cx q1[0],q0[4];
h q0[4];
sdg q0[3];
cx q1[0],q0[3];
s q0[3];
sdg q0[2];
cx q1[0],q0[2];
s q0[2];
cx q1[0],q0[1];
h q0[0];
cx q1[0],q0[0];
h q0[0];
h q1[0];
