OPENQASM 2.0;
include "qelib1.inc";
qreg q17[2];
x q17[0];
h q17[0];
x q17[1];
cx q17[0],q17[1];
h q17[0];
s q17[0];
h q17[1];
s q17[1];
cx q17[0],q17[1];
z q17[0];
h q17[0];
z q17[1];
s q17[1];
h q17[1];
s q17[1];
cx q17[0],q17[1];
sdg q17[0];
h q17[0];
h q17[1];
s q17[1];
cx q17[0],q17[1];
h q17[0];
s q17[0];
z q17[0];
h q17[0];
s q17[0];
z q17[1];
s q17[1];
h q17[1];
s q17[1];
cx q17[0],q17[1];
h q17[0];
s q17[0];
h q17[1];
s q17[1];
cx q17[0],q17[1];
z q17[0];
h q17[0];
s q17[0];
y q17[1];
sdg q17[1];
h q17[1];
cx q17[0],q17[1];
sdg q17[0];
h q17[0];
h q17[1];
s q17[1];
cx q17[0],q17[1];
h q17[0];
s q17[0];
x q17[0];
s q17[0];
h q17[0];
s q17[0];
x q17[1];
cx q17[0],q17[1];
z q17[0];
sdg q17[0];
h q17[0];
s q17[1];
h q17[1];
s q17[1];
cx q17[0],q17[1];
h q17[0];
s q17[0];
sdg q17[1];
h q17[1];
cx q17[0],q17[1];
h q17[0];
h q17[1];
s q17[1];
x q17[1];
h q17[1];
s q17[1];
cx q17[0],q17[1];
sdg q17[0];
h q17[0];
h q17[1];
s q17[1];
cx q17[0],q17[1];
y q17[0];
x q17[1];
h q17[1];
s q17[1];
cx q17[0],q17[1];
sdg q17[0];
h q17[0];
sdg q17[1];
h q17[1];
cx q17[0],q17[1];
h q17[0];
s q17[0];
y q17[0];
h q17[1];
cx q17[0],q17[1];
sdg q17[0];
h q17[0];
sdg q17[1];
h q17[1];
