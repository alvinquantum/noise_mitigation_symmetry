OPENQASM 2.0;
include "qelib1.inc";
qreg q33[2];
y q33[0];
sdg q33[0];
h q33[0];
x q33[1];
h q33[1];
s q33[1];
cx q33[0],q33[1];
sdg q33[0];
h q33[0];
sdg q33[1];
h q33[1];
cx q33[0],q33[1];
h q33[0];
s q33[0];
z q33[0];
h q33[0];
h q33[0];
h q33[1];
s q33[1];
h q33[1];
y q33[1];
s q33[1];
cx q33[0],q33[1];
h q33[0];
s q33[0];
sdg q33[1];
h q33[1];
cx q33[0],q33[1];
