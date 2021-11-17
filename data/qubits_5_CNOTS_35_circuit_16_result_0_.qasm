OPENQASM 2.0;
include "qelib1.inc";
gate circuit-9 q0,q1,q2,q3,q4 { s q0; h q0; h q1; cx q1,q2; h q2; s q3; h q3; h q4; swap q4,q0; cx q1,q0; h q0; s q0; cx q3,q4; swap q3,q0; h q0; cx q3,q2; s q2; cx q4,q1; h q1; cx q1,q4; cx q3,q1; h q1; cx q1,q3; s q1; swap q1,q0; cx q1,q2; cx q1,q0; z q1; h q1; s q2; cx q0,q2; h q2; cx q2,q0; z q0; swap q1,q0; h q2; s q2; h q2; cx q1,q2; h q2; s q3; h q3; cx q3,q1; cx q0,q1; rz(0.59838674) q3; s q3; swap q3,q0; cx q3,q0; s q0; h q0; x q4; cx q1,q4; h q4; cx q4,q1; x q1; s q4; h q4; cx q4,q3; cx q3,q2; h q2; cx q2,q3; h q2; y q3; s q3; h q3; rz(4.7515587) q4; cx q4,q0; cx q0,q2; h q2; cx q2,q0; z q0; s q0; rz(3.9213212) q0; s q2; h q2; s q2; h q4; cx q2,q4; h q4; cx q4,q2; z q2; h q2; h q4; s q4; swap q4,q0; s q0; cx q1,q4; s q1; rz(0.71250963) q4; cx q2,q4; h q2; cx q0,q2; cx q3,q4; s q3; cx q0,q3; cx q0,q1; h q1; cx q1,q0; s q1; h q1; cx q2,q1; rz(4.4089643) q2; s q2; h q2; h q3; cx q1,q3; h q3; }
qreg q0[5];
qreg q1[1];
rx(0.22930826) q0[0];
ry(1.9023143) q0[0];
rz(1.6904692) q0[0];
rx(0.21082418) q0[1];
ry(1.3605216) q0[1];
rz(1.1520803) q0[1];
rx(0.37022384) q0[2];
ry(1.1711564) q0[2];
rz(1.395833) q0[2];
rx(1.0613023) q0[3];
ry(0.053758559) q0[3];
rz(1.06075) q0[3];
rx(1.7076527) q0[4];
ry(1.2978275) q0[4];
rz(0.80851167) q0[4];
h q1[0];
barrier q0[0],q0[1],q0[2],q0[3],q0[4],q1[0];
sdg q0[3];
cx q1[0],q0[3];
s q0[3];
sdg q0[1];
cx q1[0],q0[1];
s q0[1];
barrier q0[0],q0[1],q0[2],q0[3],q0[4],q1[0];
circuit-9 q0[0],q0[1],q0[2],q0[3],q0[4];
barrier q0[0],q0[1],q0[2],q0[3],q0[4],q1[0];
cx q1[0],q0[4];
cx q1[0],q0[3];
cx q1[0],q0[2];
cx q1[0],q0[1];
cx q1[0],q0[0];
h q1[0];
