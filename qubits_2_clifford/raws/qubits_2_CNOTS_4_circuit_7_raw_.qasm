OPENQASM 2.0;
include "qelib1.inc";
qreg q6[2];
z q6[0];
x q6[1];
s q6[1];
h q6[1];
s q6[1];
cx q6[0],q6[1];
y q6[0];
s q6[0];
h q6[1];
s q6[1];
z q6[1];
s q6[1];
cx q6[0],q6[1];
sdg q6[0];
h q6[0];
x q6[0];
s q6[0];
sdg q6[1];
h q6[1];
z q6[1];
s q6[1];
cx q6[0],q6[1];
h q6[0];
s q6[0];
sdg q6[1];
h q6[1];
cx q6[0],q6[1];
