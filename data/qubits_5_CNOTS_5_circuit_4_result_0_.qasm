OPENQASM 2.0;
include "qelib1.inc";
gate circuit-7031 q0,q1,q2,q3,q4 { h q0; s q1; rz(4.7156) q1; h q3; swap q1,q4; cx q1,q4; rz(5.138882) q1; cx q0,q1; s q0; rz(5.0405047) q1; cx q2,q1; x q1; s q2; h q2; cx q3,q2; swap q4,q0; rz(1.2647672) q0; cx q0,q4; rz(3.2526346) q0; h q0; }
qreg q4[5];
qreg q5[1];
rx(0.97193365) q4[0];
ry(1.8913875) q4[0];
rz(0.93311842) q4[0];
rx(1.4517256) q4[1];
ry(0.0091718565) q4[1];
rz(0.58486132) q4[1];
rx(1.1612442) q4[2];
ry(0.22485634) q4[2];
rz(1.5110525) q4[2];
rx(1.5884188) q4[3];
ry(1.6714553) q4[3];
rz(0.71745491) q4[3];
rx(1.3221098) q4[4];
ry(0.02738419) q4[4];
rz(1.5776435) q4[4];
h q5[0];
barrier q4[0],q4[1],q4[2],q4[3],q4[4],q5[0];
sdg q4[3];
cx q5[0],q4[3];
s q4[3];
sdg q4[2];
cx q5[0],q4[2];
s q4[2];
cx q5[0],q4[0];
barrier q4[0],q4[1],q4[2],q4[3],q4[4],q5[0];
circuit-7031 q4[0],q4[1],q4[2],q4[3],q4[4];
barrier q4[0],q4[1],q4[2],q4[3],q4[4],q5[0];
h q4[4];
cx q5[0],q4[4];
h q4[4];
cx q5[0],q4[3];
sdg q4[2];
cx q5[0],q4[2];
s q4[2];
cx q5[0],q4[1];
cx q5[0],q4[0];
h q5[0];
