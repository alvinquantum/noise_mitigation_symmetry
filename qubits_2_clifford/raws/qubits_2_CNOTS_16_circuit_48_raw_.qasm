OPENQASM 2.0;
include "qelib1.inc";
qreg q47[2];
s q47[1];
h q47[1];
s q47[1];
cx q47[0],q47[1];
sdg q47[0];
h q47[0];
x q47[0];
x q47[1];
cx q47[0],q47[1];
h q47[0];
sdg q47[1];
h q47[1];
cx q47[0],q47[1];
h q47[0];
x q47[1];
h q47[1];
s q47[1];
cx q47[0],q47[1];
h q47[0];
s q47[0];
h q47[0];
s q47[0];
z q47[1];
sdg q47[1];
h q47[1];
cx q47[0],q47[1];
h q47[0];
s q47[0];
z q47[0];
s q47[0];
h q47[0];
s q47[0];
h q47[1];
s q47[1];
x q47[1];
cx q47[0],q47[1];
sdg q47[0];
h q47[0];
sdg q47[1];
h q47[1];
cx q47[0],q47[1];
h q47[0];
s q47[0];
h q47[1];
s q47[1];
cx q47[0],q47[1];
sdg q47[0];
h q47[0];
z q47[1];
sdg q47[1];
h q47[1];
cx q47[0],q47[1];
sdg q47[0];
h q47[0];
z q47[0];
s q47[0];
h q47[0];
s q47[0];
h q47[1];
s q47[1];
x q47[1];
h q47[1];
cx q47[0],q47[1];
sdg q47[0];
h q47[0];
y q47[0];
h q47[1];
s q47[1];
x q47[1];
h q47[1];
cx q47[0],q47[1];
sdg q47[0];
h q47[0];
y q47[0];
sdg q47[0];
h q47[0];
y q47[1];
sdg q47[1];
h q47[1];
cx q47[0],q47[1];
sdg q47[0];
h q47[0];
h q47[1];
s q47[1];
cx q47[0],q47[1];
x q47[0];
x q47[1];
s q47[1];
cx q47[0],q47[1];
x q47[0];
y q47[1];
s q47[1];
h q47[1];
s q47[1];
cx q47[0],q47[1];
h q47[0];
s q47[0];
x q47[0];
s q47[0];
h q47[0];
s q47[0];
sdg q47[1];
h q47[1];
s q47[1];
cx q47[0],q47[1];
sdg q47[0];
h q47[0];
sdg q47[1];
h q47[1];
