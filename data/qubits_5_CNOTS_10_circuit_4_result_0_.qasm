OPENQASM 2.0;
include "qelib1.inc";
gate circuit-3516 q0,q1,q2,q3,q4 { s q0; h q0; h q1; s q2; rz(0.10542724) q2; swap q2,q0; cx q2,q1; s q1; h q1; rz(1.8444855) q1; h q3; rz(5.1373394) q3; cx q2,q3; cx q2,q0; h q0; rz(0.73314848) q2; cx q0,q2; s q0; swap q1,q0; z q2; h q3; cx q1,q3; cx q0,q1; h q3; cx q3,q0; rz(3.0302309) q3; h q3; s q3; cx q1,q4; h q4; cx q4,q1; cx q4,q0; }
qreg q2[5];
qreg q3[1];
rx(1.3218481) q2[0];
ry(1.8847614) q2[0];
rz(1.9683396) q2[0];
rx(1.1805738) q2[1];
ry(0.015691751) q2[1];
rz(0.36351974) q2[1];
rx(1.8454065) q2[2];
ry(0.46361192) q2[2];
rz(1.5652973) q2[2];
rx(1.5005587) q2[3];
ry(1.5457862) q2[3];
rz(1.7235176) q2[3];
rx(1.3144547) q2[4];
ry(1.9587225) q2[4];
rz(0.36642626) q2[4];
h q3[0];
barrier q2[0],q2[1],q2[2],q2[3],q2[4],q3[0];
sdg q2[4];
cx q3[0],q2[4];
s q2[4];
sdg q2[0];
cx q3[0],q2[0];
s q2[0];
barrier q2[0],q2[1],q2[2],q2[3],q2[4],q3[0];
circuit-3516 q2[0],q2[1],q2[2],q2[3],q2[4];
barrier q2[0],q2[1],q2[2],q2[3],q2[4],q3[0];
cx q3[0],q2[4];
sdg q2[3];
cx q3[0],q2[3];
s q2[3];
h q2[2];
cx q3[0],q2[2];
h q2[2];
cx q3[0],q2[1];
sdg q2[0];
cx q3[0],q2[0];
s q2[0];
h q3[0];
