OPENQASM 2.0;
include "qelib1.inc";
qreg q34[2];
z q34[0];
s q34[0];
h q34[0];
s q34[0];
y q34[1];
sdg q34[1];
h q34[1];
cx q34[0],q34[1];
x q34[0];
h q34[1];
s q34[1];
cx q34[0],q34[1];
h q34[0];
s q34[0];
h q34[1];
s q34[1];
