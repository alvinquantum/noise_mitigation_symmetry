OPENQASM 2.0;
include "qelib1.inc";
gate circuit-3508 q0,q1,q2,q3,q4 { s q0; h q0; s q1; rz(2.1864052) q1; h q2; rz(2.9669408) q2; cx q1,q2; h q2; h q3; rz(0.72696354) q3; swap q3,q0; rz(1.7203774) q4; h q4; cx q1,q4; cx q1,q0; h q0; rz(2.393222) q0; s q0; cx q3,q1; h q1; cx q1,q3; h q4; }
qreg q2[5];
qreg q3[1];
rx(0.33741707) q2[0];
ry(0.41092977) q2[0];
rz(0.8106614) q2[0];
rx(1.0627708) q2[1];
ry(1.7245081) q2[1];
rz(1.689058) q2[1];
rx(1.4022236) q2[2];
ry(0.52091517) q2[2];
rz(0.058633987) q2[2];
rx(1.6240334) q2[3];
ry(0.23917698) q2[3];
rz(0.70757976) q2[3];
rx(0.24683293) q2[4];
ry(1.9004199) q2[4];
rz(1.151599) q2[4];
h q3[0];
barrier q2[0],q2[1],q2[2],q2[3],q2[4],q3[0];
h q2[4];
cx q3[0],q2[4];
h q2[4];
cx q3[0],q2[2];
h q2[1];
cx q3[0],q2[1];
h q2[1];
h q2[0];
cx q3[0],q2[0];
h q2[0];
barrier q2[0],q2[1],q2[2],q2[3],q2[4],q3[0];
circuit-3508 q2[0],q2[1],q2[2],q2[3],q2[4];
barrier q2[0],q2[1],q2[2],q2[3],q2[4],q3[0];
h q2[4];
cx q3[0],q2[4];
h q2[4];
cx q3[0],q2[3];
cx q3[0],q2[2];
h q2[1];
cx q3[0],q2[1];
h q2[1];
h q3[0];
