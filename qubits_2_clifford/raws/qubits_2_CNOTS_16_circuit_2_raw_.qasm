OPENQASM 2.0;
include "qelib1.inc";
qreg q1[2];
y q1[0];
s q1[0];
h q1[0];
s q1[0];
y q1[1];
cx q1[0],q1[1];
y q1[0];
sdg q1[0];
h q1[0];
z q1[1];
s q1[1];
h q1[1];
s q1[1];
cx q1[0],q1[1];
h q1[0];
s q1[0];
sdg q1[1];
h q1[1];
cx q1[0],q1[1];
z q1[0];
sdg q1[0];
h q1[0];
h q1[1];
cx q1[0],q1[1];
h q1[0];
s q1[0];
sdg q1[1];
h q1[1];
cx q1[0],q1[1];
x q1[1];
h q1[1];
cx q1[0],q1[1];
sdg q1[0];
h q1[0];
sdg q1[0];
h q1[0];
x q1[1];
h q1[1];
cx q1[0],q1[1];
sdg q1[0];
h q1[0];
h q1[1];
s q1[1];
cx q1[0],q1[1];
y q1[0];
s q1[0];
h q1[0];
s q1[0];
y q1[1];
s q1[1];
cx q1[0],q1[1];
sdg q1[0];
h q1[0];
h q1[1];
s q1[1];
cx q1[0],q1[1];
h q1[0];
s q1[0];
z q1[0];
sdg q1[0];
h q1[0];
z q1[1];
s q1[1];
h q1[1];
s q1[1];
cx q1[0],q1[1];
sdg q1[0];
h q1[0];
sdg q1[1];
h q1[1];
cx q1[0],q1[1];
h q1[0];
s q1[0];
x q1[0];
s q1[1];
cx q1[0],q1[1];
x q1[0];
s q1[0];
h q1[0];
s q1[0];
x q1[1];
h q1[1];
s q1[1];
cx q1[0],q1[1];
sdg q1[0];
h q1[0];
sdg q1[1];
h q1[1];
cx q1[0],q1[1];
z q1[0];
h q1[1];
s q1[1];
z q1[1];
h q1[1];
s q1[1];
cx q1[0],q1[1];
sdg q1[0];
h q1[0];
sdg q1[1];
h q1[1];
