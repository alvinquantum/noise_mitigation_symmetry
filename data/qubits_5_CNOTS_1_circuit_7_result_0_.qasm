OPENQASM 2.0;
include "qelib1.inc";
gate circuit-9 q0,q1,q2,q3,q4 { rz(1.4058162) q0; s q0; h q0; s q1; s q2; s q3; h q3; rz(1.9193987) q3; s q4; rz(5.9619295) q4; h q4; rz(2.2302651) q4; cx q4,q0; rz(4.1219015) q0; }
qreg q0[5];
qreg q1[1];
rx(1.7985771) q0[0];
ry(1.7063185) q0[0];
rz(0.93021459) q0[0];
rx(1.4594421) q0[1];
ry(1.2540488) q0[1];
rz(0.74404998) q0[1];
rx(0.052989781) q0[2];
ry(1.2345814) q0[2];
rz(1.5247219) q0[2];
rx(0.71497103) q0[3];
ry(1.9350274) q0[3];
rz(1.6476246) q0[3];
rx(0.9275643) q0[4];
ry(1.5232918) q0[4];
rz(1.5215196) q0[4];
h q1[0];
barrier q0[0],q0[1],q0[2],q0[3],q0[4],q1[0];
sdg q0[3];
cx q1[0],q0[3];
s q0[3];
cx q1[0],q0[2];
sdg q0[1];
cx q1[0],q0[1];
s q0[1];
barrier q0[0],q0[1],q0[2],q0[3],q0[4],q1[0];
circuit-9 q0[0],q0[1],q0[2],q0[3],q0[4];
barrier q0[0],q0[1],q0[2],q0[3],q0[4],q1[0];
h q0[3];
cx q1[0],q0[3];
h q0[3];
sdg q0[2];
cx q1[0],q0[2];
s q0[2];
cx q1[0],q0[1];
h q1[0];
