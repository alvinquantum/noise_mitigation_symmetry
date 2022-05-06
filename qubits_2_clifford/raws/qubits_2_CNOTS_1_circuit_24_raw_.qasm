OPENQASM 2.0;
include "qelib1.inc";
qreg q23[2];
y q23[0];
s q23[0];
s q23[1];
cx q23[0],q23[1];
