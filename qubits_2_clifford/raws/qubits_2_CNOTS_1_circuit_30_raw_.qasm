OPENQASM 2.0;
include "qelib1.inc";
qreg q29[2];
x q29[0];
sdg q29[0];
h q29[0];
h q29[1];
cx q29[0],q29[1];
h q29[0];
s q29[0];
sdg q29[1];
h q29[1];
