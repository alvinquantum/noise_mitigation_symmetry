OPENQASM 2.0;
include "qelib1.inc";
qreg q32[2];
z q32[0];
h q32[0];
s q32[0];
h q32[1];
cx q32[0],q32[1];
h q32[0];
s q32[0];
h q32[1];
s q32[1];
cx q32[0],q32[1];
