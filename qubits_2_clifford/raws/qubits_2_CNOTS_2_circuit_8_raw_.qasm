OPENQASM 2.0;
include "qelib1.inc";
qreg q7[2];
y q7[0];
h q7[0];
y q7[1];
s q7[1];
cx q7[0],q7[1];
x q7[0];
s q7[0];
h q7[0];
s q7[0];
h q7[1];
s q7[1];
cx q7[0],q7[1];
sdg q7[0];
h q7[0];
sdg q7[1];
h q7[1];
