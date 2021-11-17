OPENQASM 2.0;
include "qelib1.inc";
gate circuit-9 q0,q1,q2,q3,q4 { s q0; h q0; s q1; h q2; cx q0,q2; h q2; h q3; cx q0,q3; cx q2,q3; rz(4.6311665) q3; h q3; s q4; h q4; cx q4,q0; rz(2.594758) q0; cx q0,q1; h q1; cx q1,q0; s q1; rz(5.8417604) q1; h q1; s q1; cx q4,q1; cx q1,q2; rz(4.9748532) q1; h q1; s q1; cx q2,q4; x q2; swap q2,q0; cx q2,q0; h q0; s q4; swap q4,q1; cx q4,q1; cx q3,q4; s q3; h q3; swap q3,q1; cx q1,q3; h q1; cx q1,q2; swap q1,q0; rz(1.8858806) q3; y q3; h q3; cx q3,q2; cx q2,q4; h q3; h q4; cx q4,q2; cx q1,q4; x q2; cx q3,q1; cx q1,q0; h q0; cx q0,q1; y q1; s q1; h q1; h q3; cx q0,q3; h q4; cx q4,q0; z q0; h q0; s q0; swap q4,q3; cx q4,q3; z q3; s q3; y q4; cx q0,q4; cx q1,q0; h q1; cx q2,q0; cx q0,q3; cx q1,q2; s q2; h q3; h q4; cx q1,q4; }
qreg q0[5];
qreg q1[1];
rx(0.50991779) q0[0];
ry(0.79154446) q0[0];
rz(0.29263973) q0[0];
rx(0.68377235) q0[1];
ry(1.5753132) q0[1];
rz(1.9587539) q0[1];
rx(0.52252015) q0[2];
ry(0.10665525) q0[2];
rz(1.306534) q0[2];
rx(0.71415333) q0[3];
ry(1.0243621) q0[3];
rz(1.5312986) q0[3];
rx(0.34606009) q0[4];
ry(0.2111607) q0[4];
rz(1.1168149) q0[4];
h q1[0];
barrier q0[0],q0[1],q0[2],q0[3],q0[4],q1[0];
cx q1[0],q0[4];
cx q1[0],q0[3];
sdg q0[2];
cx q1[0],q0[2];
s q0[2];
cx q1[0],q0[0];
barrier q0[0],q0[1],q0[2],q0[3],q0[4],q1[0];
circuit-9 q0[0],q0[1],q0[2],q0[3],q0[4];
barrier q0[0],q0[1],q0[2],q0[3],q0[4],q1[0];
sdg q0[4];
cx q1[0],q0[4];
s q0[4];
h q0[3];
cx q1[0],q0[3];
h q0[3];
cx q1[0],q0[2];
cx q1[0],q0[1];
h q0[0];
cx q1[0],q0[0];
h q0[0];
h q1[0];
