OPENQASM 2.0;
include "qelib1.inc";
qreg q10[2];
x q10[0];
z q10[1];
s q10[1];
h q10[1];
s q10[1];
cx q10[0],q10[1];
h q10[0];
s q10[0];
sdg q10[1];
h q10[1];
cx q10[0],q10[1];
h q10[0];
z q10[1];
h q10[1];
s q10[1];
cx q10[0],q10[1];
z q10[0];
h q10[0];
s q10[0];
y q10[0];
h q10[0];
sdg q10[1];
h q10[1];
x q10[1];
sdg q10[1];
h q10[1];
x q10[1];
s q10[1];
x q10[1];
s q10[1];
cx q10[0],q10[1];
h q10[0];
s q10[0];
z q10[0];
s q10[1];
h q10[1];
s q10[1];
cx q10[0],q10[1];
sdg q10[0];
h q10[0];
sdg q10[1];
h q10[1];
cx q10[0],q10[1];
x q10[0];
sdg q10[0];
h q10[0];
h q10[1];
s q10[1];
x q10[1];
h q10[1];
cx q10[0],q10[1];
h q10[0];
s q10[0];
h q10[1];
s q10[1];
cx q10[0],q10[1];
y q10[0];
h q10[0];
z q10[1];
s q10[1];
h q10[1];
s q10[1];
cx q10[0],q10[1];
sdg q10[0];
h q10[0];
s q10[0];
s q10[1];
cx q10[0],q10[1];
h q10[0];
s q10[0];
y q10[0];
sdg q10[0];
h q10[0];
z q10[1];
h q10[1];
cx q10[0],q10[1];
sdg q10[0];
h q10[0];
h q10[1];
s q10[1];
cx q10[0],q10[1];
x q10[0];
s q10[0];
h q10[0];
s q10[0];
sdg q10[1];
h q10[1];
cx q10[0],q10[1];
sdg q10[0];
h q10[0];
x q10[0];
s q10[0];
h q10[1];
s q10[1];
cx q10[0],q10[1];
sdg q10[0];
h q10[0];
y q10[0];
h q10[0];
sdg q10[1];
h q10[1];
z q10[1];
s q10[1];
h q10[1];
s q10[1];
cx q10[0],q10[1];
h q10[0];
s q10[0];
x q10[0];
h q10[0];
s q10[0];
sdg q10[1];
h q10[1];
y q10[1];
s q10[1];
cx q10[0],q10[1];
h q10[0];
s q10[0];
sdg q10[1];
h q10[1];
