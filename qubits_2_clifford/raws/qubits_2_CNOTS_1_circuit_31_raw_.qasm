OPENQASM 2.0;
include "qelib1.inc";
qreg q30[2];
h q30[0];
x q30[1];
cx q30[0],q30[1];
