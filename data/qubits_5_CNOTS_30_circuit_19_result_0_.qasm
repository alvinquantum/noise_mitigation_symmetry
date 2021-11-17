OPENQASM 2.0;
include "qelib1.inc";
gate circuit-14049 q0,q1,q2,q3,q4 { s q1; h q2; cx q0,q2; s q3; cx q3,q1; rz(4.367391) q4; cx q4,q0; cx q1,q0; cx q0,q3; h q0; rz(3.3748299) q1; s q1; h q1; s q1; h q3; cx q3,q1; s q4; cx q2,q4; cx q1,q4; h q1; h q2; s q2; cx q2,q1; cx q4,q3; s q3; h q3; s q3; cx q1,q3; cx q3,q2; cx q1,q2; h q2; cx q2,q1; y q1; s q2; h q2; s q2; x q2; h q2; cx q2,q1; rz(3.7002737) q1; h q1; rz(4.0767778) q3; y q3; cx q2,q3; h q3; cx q3,q2; h q3; h q4; cx q4,q0; cx q0,q2; h q0; s q0; cx q2,q4; rz(1.1371286) q2; x q2; s q2; h q2; h q4; cx q4,q0; cx q0,q1; h q1; cx q1,q0; s q0; h q1; s q1; swap q4,q1; cx q4,q1; h q1; cx q1,q4; s q1; h q1; cx q1,q3; x q1; h q3; x q3; h q3; cx q3,q2; z q4; s q4; cx q0,q4; cx q0,q1; h q1; cx q1,q0; s q4; h q4; s q4; }
qreg q8[5];
qreg q9[1];
rx(1.8494546) q8[0];
ry(1.7310503) q8[0];
rz(1.7218763) q8[0];
rx(1.7632946) q8[1];
ry(1.0945733) q8[1];
rz(1.6296087) q8[1];
rx(0.38923057) q8[2];
ry(1.20505) q8[2];
rz(1.8013722) q8[2];
rx(1.5360003) q8[3];
ry(0.57473384) q8[3];
rz(1.8898845) q8[3];
rx(1.189399) q8[4];
ry(1.3885519) q8[4];
rz(1.9755341) q8[4];
h q9[0];
barrier q8[0],q8[1],q8[2],q8[3],q8[4],q9[0];
cx q9[0],q8[3];
sdg q8[1];
cx q9[0],q8[1];
s q8[1];
sdg q8[0];
cx q9[0],q8[0];
s q8[0];
barrier q8[0],q8[1],q8[2],q8[3],q8[4],q9[0];
circuit-14049 q8[0],q8[1],q8[2],q8[3],q8[4];
barrier q8[0],q8[1],q8[2],q8[3],q8[4],q9[0];
cx q9[0],q8[4];
h q8[3];
cx q9[0],q8[3];
h q8[3];
cx q9[0],q8[2];
sdg q8[1];
cx q9[0],q8[1];
s q8[1];
sdg q8[0];
cx q9[0],q8[0];
s q8[0];
h q9[0];
