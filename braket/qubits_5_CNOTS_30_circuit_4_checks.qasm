OPENQASM 2.0;
include "qelib1.inc";
qreg q[6];
creg c[1];
rx(1.8053908) q[0];
ry(1.5731695) q[0];
rz(0.95543868) q[0];
rx(1.8592313) q[1];
ry(1.75681) q[1];
rz(1.7498729) q[1];
rx(1.2571594) q[2];
ry(0.049380882) q[2];
rz(1.5357697) q[2];
rx(0.10517681) q[3];
ry(0.16478934) q[3];
rz(1.2837496) q[3];
rx(1.7225831) q[4];
ry(0.9939463) q[4];
rz(0.87938714) q[4];
barrier q[0],q[1],q[2],q[3],q[4],q[5];
h q[5];
barrier q[0],q[1],q[2],q[3],q[4],q[5];
sdg q[0];
cx q[5],q[3];
cx q[5],q[2];
cx q[5],q[0];
s q[0];
barrier q[0],q[1],q[2],q[3],q[4],q[5];
h q[0];
h q[1];
s q[2];
s q[3];
cx q[3],q[1];
h q[1];
s q[1];
cx q[3],q[0];
h q[0];
cx q[1],q[0];
cx q[0],q[2];
cx q[2],q[1];
h q[1];
z q[2];
rz(4.1878145) q[2];
s q[2];
rz(0.18837346) q[4];
h q[4];
cx q[3],q[4];
h q[4];
cx q[4],q[3];
z q[3];
h q[3];
s q[3];
h q[4];
s q[4];
cx q[4],q[1];
cx q[1],q[0];
cx q[0],q[4];
z q[0];
h q[1];
s q[1];
cx q[4],q[1];
s q[1];
h q[1];
s q[1];
x q[1];
s q[1];
h q[1];
cx q[3],q[1];
y q[4];
cx q[0],q[4];
h q[0];
cx q[4],q[2];
h q[2];
cx q[2],q[4];
cx q[1],q[4];
h q[2];
cx q[2],q[0];
h q[0];
cx q[0],q[2];
cx q[4],q[3];
rz(2.9714438) q[3];
s q[3];
h q[3];
cx q[3],q[1];
cx q[1],q[2];
rz(0.95384789) q[2];
cx q[2],q[3];
z q[2];
s q[2];
h q[2];
s q[2];
s q[3];
cx q[3],q[0];
h q[0];
s q[0];
cx q[3],q[1];
h q[1];
cx q[1],q[3];
cx q[0],q[1];
h q[0];
s q[0];
s q[0];
h q[0];
s q[1];
h q[1];
s q[1];
rz(0.65267227) q[1];
s q[3];
h q[3];
cx q[3],q[0];
s q[3];
y q[4];
cx q[4],q[2];
cx q[2],q[0];
cx q[0],q[4];
barrier q[0],q[1],q[2],q[3],q[4],q[5];
sdg q[0];
h q[1];
h q[3];
cx q[5],q[4];
cx q[5],q[3];
cx q[5],q[2];
cx q[5],q[1];
cx q[5],q[0];
s q[0];
h q[1];
h q[3];
barrier q[0],q[1],q[2],q[3],q[4],q[5];
h q[5];
barrier q[0],q[1],q[2],q[3],q[4],q[5];
measure q[5] -> c[0];
