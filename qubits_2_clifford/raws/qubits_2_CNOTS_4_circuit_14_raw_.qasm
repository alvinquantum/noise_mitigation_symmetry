OPENQASM 2.0;
include "qelib1.inc";
qreg q13[2];
z q13[0];
y q13[1];
sdg q13[1];
h q13[1];
cx q13[0],q13[1];
sdg q13[0];
h q13[0];
sdg q13[1];
h q13[1];
cx q13[0],q13[1];
h q13[0];
s q13[0];
h q13[1];
s q13[1];
cx q13[0],q13[1];
s q13[0];
y q13[1];
cx q13[0],q13[1];
sdg q13[0];
h q13[0];
