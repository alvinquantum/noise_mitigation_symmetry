OPENQASM 2.0;
include "qelib1.inc";
qreg q24[2];
z q24[0];
h q24[0];
x q24[1];
sdg q24[1];
h q24[1];
cx q24[0],q24[1];
sdg q24[0];
h q24[0];
h q24[1];
s q24[1];
cx q24[0],q24[1];
