OPENQASM 2.0;
include "qelib1.inc";
qreg q45[2];
y q45[0];
x q45[1];
sdg q45[1];
h q45[1];
cx q45[0],q45[1];
sdg q45[0];
h q45[0];
