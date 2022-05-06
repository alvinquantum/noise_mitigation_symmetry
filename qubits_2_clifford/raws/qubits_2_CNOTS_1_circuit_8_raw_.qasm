OPENQASM 2.0;
include "qelib1.inc";
qreg q7[2];
z q7[0];
s q7[0];
cx q7[0],q7[1];
h q7[0];
s q7[0];
