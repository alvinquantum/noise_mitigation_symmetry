OPENQASM 2.0;
include "qelib1.inc";
gate circuit-3520 q0,q1,q2,q3,q4 { s q2; h q2; s q2; swap q1,q2; h q2; cx q1,q3; cx q0,q1; s q0; s q3; cx q3,q0; s q4; h q4; cx q4,q1; y q1; cx q4,q2; h q2; cx q2,q4; cx q0,q4; h q0; s q0; swap q2,q0; s q0; h q0; cx q4,q3; h q3; cx q3,q2; cx q0,q3; x q2; s q2; h q2; s q2; h q3; cx q3,q0; y q0; h q0; s q0; cx q1,q0; s q1; h q1; s q3; h q3; z q3; cx q3,q2; rz(5.2870761) q3; y q4; s q4; cx q0,q4; h q4; cx q4,q0; cx q2,q0; cx q0,q3; z q0; h q0; s q0; s q2; cx q2,q1; s q4; h q4; cx q3,q4; h q3; rz(2.7515946) q4; cx q1,q4; s q1; swap q3,q1; cx q1,q3; h q1; cx q4,q2; s q2; cx q3,q2; h q2; cx q2,q3; s q2; h q2; swap q2,q1; cx q1,q2; h q1; y q1; h q1; cx q0,q1; h q1; cx q1,q0; s q1; h q1; s q1; y q2; rz(1.5349647) q3; x q3; x q4; h q4; s q4; rz(5.2799934) q4; cx q4,q3; cx q3,q0; cx q0,q4; s q0; s q3; h q3; s q4; swap q4,q1; cx q4,q1; s q1; h q1; cx q3,q4; h q3; swap q3,q1; cx q4,q2; rz(5.9238976) q2; h q2; cx q2,q4; h q2; cx q2,q1; s q1; h q1; cx q3,q2; h q2; cx q2,q3; h q2; cx q1,q2; y q1; swap q1,q0; h q2; s q2; s q2; x q3; s q3; cx q3,q2; y q4; cx q1,q4; }
qreg q2[5];
qreg q3[1];
rx(0.2048235) q2[0];
ry(0.40065873) q2[0];
rz(1.8836264) q2[0];
rx(1.5917464) q2[1];
ry(1.8953662) q2[1];
rz(1.1235679) q2[1];
rx(1.2386094) q2[2];
ry(1.372171) q2[2];
rz(1.1023725) q2[2];
rx(1.9231908) q2[3];
ry(0.43221098) q2[3];
rz(0.8690936) q2[3];
rx(1.2923752) q2[4];
ry(0.98716854) q2[4];
rz(0.028885361) q2[4];
h q3[0];
barrier q2[0],q2[1],q2[2],q2[3],q2[4],q3[0];
cx q3[0],q2[4];
cx q3[0],q2[3];
cx q3[0],q2[2];
h q2[1];
cx q3[0],q2[1];
h q2[1];
h q2[0];
cx q3[0],q2[0];
h q2[0];
barrier q2[0],q2[1],q2[2],q2[3],q2[4],q3[0];
circuit-3520 q2[0],q2[1],q2[2],q2[3],q2[4];
barrier q2[0],q2[1],q2[2],q2[3],q2[4],q3[0];
sdg q2[4];
cx q3[0],q2[4];
s q2[4];
sdg q2[3];
cx q3[0],q2[3];
s q2[3];
h q2[2];
cx q3[0],q2[2];
h q2[2];
sdg q2[1];
cx q3[0],q2[1];
s q2[1];
sdg q2[0];
cx q3[0],q2[0];
s q2[0];
h q3[0];
