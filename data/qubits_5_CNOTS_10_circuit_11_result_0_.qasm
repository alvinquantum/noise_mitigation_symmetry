OPENQASM 2.0;
include "qelib1.inc";
gate circuit-14037 q0,q1,q2,q3,q4 { rz(0.77857653) q0; s q0; s q1; h q1; rz(5.4173787) q1; swap q3,q0; cx q3,q0; cx q1,q3; rz(1.213539) q3; cx q4,q3; cx q3,q2; h q2; cx q2,q3; s q2; h q2; s q2; cx q2,q1; rz(4.7121954) q1; y q3; rz(5.2248217) q3; s q4; h q4; cx q4,q0; cx q1,q0; cx q0,q2; y q0; s q1; h q2; cx q2,q1; h q4; }
qreg q8[5];
qreg q9[1];
rx(0.029000793) q8[0];
ry(0.055995761) q8[0];
rz(0.2139499) q8[0];
rx(0.084264055) q8[1];
ry(0.16712797) q8[1];
rz(1.9742417) q8[1];
rx(0.2975703) q8[2];
ry(0.52605725) q8[2];
rz(1.4112286) q8[2];
rx(0.37164092) q8[3];
ry(0.52331831) q8[3];
rz(0.64778265) q8[3];
rx(0.28663507) q8[4];
ry(0.30411107) q8[4];
rz(1.8634953) q8[4];
h q9[0];
barrier q8[0],q8[1],q8[2],q8[3],q8[4],q9[0];
sdg q8[3];
cx q9[0],q8[3];
s q8[3];
cx q9[0],q8[2];
sdg q8[1];
cx q9[0],q8[1];
s q8[1];
barrier q8[0],q8[1],q8[2],q8[3],q8[4],q9[0];
circuit-14037 q8[0],q8[1],q8[2],q8[3],q8[4];
barrier q8[0],q8[1],q8[2],q8[3],q8[4],q9[0];
sdg q8[4];
cx q9[0],q8[4];
s q8[4];
h q8[3];
cx q9[0],q8[3];
h q8[3];
h q8[2];
cx q9[0],q8[2];
h q8[2];
h q8[1];
cx q9[0],q8[1];
h q8[1];
h q8[0];
cx q9[0],q8[0];
h q8[0];
h q9[0];
