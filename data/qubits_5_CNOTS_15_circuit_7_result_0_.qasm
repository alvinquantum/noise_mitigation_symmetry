OPENQASM 2.0;
include "qelib1.inc";
gate circuit-14029 q0,q1,q2,q3,q4 { s q0; h q0; s q0; h q1; cx q0,q1; s q1; h q1; rz(0.4261564) q2; h q2; s q4; cx q0,q4; cx q3,q0; cx q0,q2; h q2; cx q2,q0; y q0; s q0; s q2; h q2; h q3; cx q1,q3; cx q2,q1; s q2; rz(0.35235182) q2; h q2; rz(1.6093964) q2; z q3; s q3; cx q0,q3; h q3; cx q4,q1; z q1; cx q1,q0; s q1; h q1; s q1; s q4; swap q4,q2; cx q4,q2; h q2; cx q2,q4; h q2; z q2; s q2; y q4; s q4; cx q4,q2; cx q2,q0; rz(1.7033556) q4; cx q0,q4; rz(3.1375272) q0; x q0; h q4; }
qreg q8[5];
qreg q9[1];
rx(1.8359261) q8[0];
ry(0.17912174) q8[0];
rz(1.3853062) q8[0];
rx(1.7293317) q8[1];
ry(0.77397523) q8[1];
rz(1.7432701) q8[1];
rx(0.56297382) q8[2];
ry(1.7623806) q8[2];
rz(1.3043108) q8[2];
rx(0.012226651) q8[3];
ry(3/(14*pi)) q8[3];
rz(0.28610308) q8[3];
rx(0.73724499) q8[4];
ry(0.18722127) q8[4];
rz(1.3862961) q8[4];
h q9[0];
barrier q8[0],q8[1],q8[2],q8[3],q8[4],q9[0];
h q8[4];
cx q9[0],q8[4];
h q8[4];
sdg q8[3];
cx q9[0],q8[3];
s q8[3];
h q8[2];
cx q9[0],q8[2];
h q8[2];
cx q9[0],q8[1];
h q8[0];
cx q9[0],q8[0];
h q8[0];
barrier q8[0],q8[1],q8[2],q8[3],q8[4],q9[0];
circuit-14029 q8[0],q8[1],q8[2],q8[3],q8[4];
barrier q8[0],q8[1],q8[2],q8[3],q8[4],q9[0];
cx q9[0],q8[4];
cx q9[0],q8[3];
sdg q8[2];
cx q9[0],q8[2];
s q8[2];
cx q9[0],q8[1];
h q8[0];
cx q9[0],q8[0];
h q8[0];
h q9[0];
