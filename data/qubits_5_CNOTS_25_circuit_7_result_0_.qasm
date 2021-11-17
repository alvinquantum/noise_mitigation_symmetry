OPENQASM 2.0;
include "qelib1.inc";
gate circuit-3520 q0,q1,q2,q3,q4 { h q0; s q0; h q1; h q2; swap q2,q0; rz(3.3107102) q2; cx q2,q1; cx q2,q0; s q0; cx q0,q1; rz(3.4235185) q1; h q1; cx q1,q0; h q1; rz(0.85071713) q1; h q3; s q3; s q4; h q4; cx q4,q3; cx q3,q2; cx q2,q4; h q2; s q3; h q3; rz(0.52053748) q4; h q4; s q4; cx q4,q3; cx q3,q0; cx q0,q4; z q0; h q0; h q3; cx q3,q1; s q4; h q4; s q4; cx q1,q4; h q1; cx q4,q3; h q3; cx q3,q1; y q1; swap q1,q0; cx q1,q2; cx q1,q0; h q0; s q0; h q3; z q3; h q3; h q4; cx q4,q1; cx q1,q3; h q3; cx q3,q1; x q1; h q4; swap q4,q0; cx q2,q4; cx q0,q4; s q0; h q0; s q2; h q2; s q2; cx q4,q3; h q3; cx q3,q4; cx q3,q2; cx q2,q0; rz(5.2914631) q3; }
qreg q2[5];
qreg q3[1];
rx(1.2766016) q2[0];
ry(1.1418008) q2[0];
rz(0.034705075) q2[0];
rx(1.8955018) q2[1];
ry(0.92504666) q2[1];
rz(0.017949658) q2[1];
rx(1.8039874) q2[2];
ry(0.67189216) q2[2];
rz(1.4130256) q2[2];
rx(1.9521317) q2[3];
ry(1.5182692) q2[3];
rz(1.7380963) q2[3];
rx(0.87710498) q2[4];
ry(1.8858958) q2[4];
rz(0.47554417) q2[4];
h q3[0];
barrier q2[0],q2[1],q2[2],q2[3],q2[4],q3[0];
cx q3[0],q2[4];
h q2[2];
cx q3[0],q2[2];
h q2[2];
sdg q2[1];
cx q3[0],q2[1];
s q2[1];
barrier q2[0],q2[1],q2[2],q2[3],q2[4],q3[0];
circuit-3520 q2[0],q2[1],q2[2],q2[3],q2[4];
barrier q2[0],q2[1],q2[2],q2[3],q2[4],q3[0];
cx q3[0],q2[4];
h q2[3];
cx q3[0],q2[3];
h q2[3];
sdg q2[2];
cx q3[0],q2[2];
s q2[2];
cx q3[0],q2[1];
cx q3[0],q2[0];
h q3[0];
