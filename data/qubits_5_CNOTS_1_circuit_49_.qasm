OPENQASM 2.0;
include "qelib1.inc";
qreg q198[5];
qreg q199[1];
h q198[4];
cx q199[0],q198[4];
h q198[4];
cx q199[0],q198[2];
h q198[1];
cx q199[0],q198[1];
h q198[1];
h q198[0];
cx q199[0],q198[0];
h q198[0];
barrier q198[0],q198[1],q198[2],q198[3],q198[4],q199[0];
s q198[0];
rz(5.733335) q198[1];
s q198[1];
s q198[2];
rz(1.6913565) q198[4];
s q198[4];
rz(1.8124867) q198[0];
swap q198[1],q198[0];
h q198[2];
rz(0.34197747) q198[3];
cx q198[3],q198[1];
s q198[2];
rz(5.2699807) q198[1];
barrier q198[0],q198[1],q198[2],q198[3],q198[4],q199[0];
h q198[4];
cx q199[0],q198[4];
h q198[4];
h q198[3];
cx q199[0],q198[3];
h q198[3];
cx q199[0],q198[2];
h q198[1];
cx q199[0],q198[1];
h q198[1];
h q198[0];
cx q199[0],q198[0];
h q198[0];
