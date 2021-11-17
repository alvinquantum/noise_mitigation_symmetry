OPENQASM 2.0;
include "qelib1.inc";
gate circuit-10518 q0,q1,q2,q3,q4 { s q1; cx q0,q2; h q2; cx q2,q0; rz(5.8832093) q0; cx q1,q2; h q2; s q3; s q4; rz(4.3327194) q4; h q4; cx q4,q3; cx q3,q0; cx q0,q4; h q0; s q0; cx q1,q4; cx q2,q4; z q2; s q2; h q2; h q3; cx q1,q3; h q3; cx q3,q1; y q1; s q3; h q3; s q4; cx q3,q4; x q3; s q3; cx q3,q1; s q4; h q4; s q4; h q4; rz(1.8442052) q4; cx q0,q4; cx q2,q0; cx q1,q0; cx q0,q3; y q0; h q0; s q0; h q1; s q1; cx q1,q3; s q2; h q3; cx q3,q1; h q3; s q4; h q4; cx q4,q2; cx q2,q1; cx q1,q4; x q1; s q2; h q2; s q2; cx q3,q2; h q3; s q3; s q4; h q4; rz(3.0892718) q4; cx q4,q2; y q2; s q2; h q2; swap q2,q0; cx q3,q4; z q3; s q3; cx q1,q3; h q4; s q4; y q4; s q4; cx q2,q4; cx q0,q2; swap q0,q3; cx q2,q1; h q1; cx q1,q2; s q1; cx q0,q1; h q1; cx q1,q0; x q0; s q1; h q1; s q1; h q3; cx q3,q1; h q4; cx q1,q4; rz(5.003716) q1; h q1; cx q4,q3; h q3; cx q1,q3; y q1; h q3; s q3; z q3; }
qreg q6[5];
qreg q7[1];
rx(1.8473397) q6[0];
ry(0.8014337) q6[0];
rz(0.7678033) q6[0];
rx(0.79188408) q6[1];
ry(0.5625485) q6[1];
rz(1.8158099) q6[1];
rx(0.76652631) q6[2];
ry(0.32734868) q6[2];
rz(0.50615705) q6[2];
rx(1.2381782) q6[3];
ry(1.7391212) q6[3];
rz(1.6902483) q6[3];
rx(1.3423658) q6[4];
ry(0.79873644) q6[4];
rz(0.95743559) q6[4];
h q7[0];
barrier q6[0],q6[1],q6[2],q6[3],q6[4],q7[0];
h q6[4];
cx q7[0],q6[4];
h q6[4];
cx q7[0],q6[3];
h q6[2];
cx q7[0],q6[2];
h q6[2];
sdg q6[1];
cx q7[0],q6[1];
s q6[1];
sdg q6[0];
cx q7[0],q6[0];
s q6[0];
barrier q6[0],q6[1],q6[2],q6[3],q6[4],q7[0];
circuit-10518 q6[0],q6[1],q6[2],q6[3],q6[4];
barrier q6[0],q6[1],q6[2],q6[3],q6[4],q7[0];
h q6[4];
cx q7[0],q6[4];
h q6[4];
cx q7[0],q6[3];
cx q7[0],q6[2];
h q6[1];
cx q7[0],q6[1];
h q6[1];
cx q7[0],q6[0];
h q7[0];
