OPENQASM 2.0;
include "qelib1.inc";
qreg q34[2];
x q34[0];
s q34[0];
sdg q34[1];
h q34[1];
cx q34[0],q34[1];
sdg q34[0];
h q34[0];
sdg q34[1];
h q34[1];
