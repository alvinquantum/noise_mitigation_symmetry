OPENQASM 2.0;
include "qelib1.inc";
qreg q22[2];
z q22[0];
sdg q22[0];
h q22[0];
h q22[1];
cx q22[0],q22[1];
