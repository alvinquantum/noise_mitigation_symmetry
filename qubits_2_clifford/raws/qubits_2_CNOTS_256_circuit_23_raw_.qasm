OPENQASM 2.0;
include "qelib1.inc";
qreg q22[2];
s q22[0];
y q22[1];
sdg q22[1];
h q22[1];
cx q22[0],q22[1];
sdg q22[0];
h q22[0];
sdg q22[1];
h q22[1];
cx q22[0],q22[1];
s q22[0];
x q22[1];
h q22[1];
cx q22[0],q22[1];
sdg q22[0];
h q22[0];
y q22[0];
h q22[1];
s q22[1];
z q22[1];
h q22[1];
s q22[1];
cx q22[0],q22[1];
sdg q22[0];
h q22[0];
z q22[0];
sdg q22[0];
h q22[0];
h q22[1];
s q22[1];
h q22[1];
s q22[1];
cx q22[0],q22[1];
h q22[0];
s q22[0];
y q22[0];
s q22[0];
h q22[0];
s q22[0];
sdg q22[1];
h q22[1];
x q22[1];
s q22[1];
cx q22[0],q22[1];
sdg q22[0];
h q22[0];
sdg q22[1];
h q22[1];
cx q22[0],q22[1];
h q22[0];
s q22[0];
h q22[0];
s q22[0];
x q22[1];
sdg q22[1];
h q22[1];
cx q22[0],q22[1];
h q22[0];
s q22[0];
sdg q22[1];
h q22[1];
cx q22[0],q22[1];
x q22[0];
s q22[0];
h q22[0];
s q22[0];
h q22[1];
s q22[1];
sdg q22[1];
h q22[1];
cx q22[0],q22[1];
sdg q22[0];
h q22[0];
sdg q22[1];
h q22[1];
cx q22[0],q22[1];
y q22[0];
h q22[0];
s q22[0];
s q22[1];
cx q22[0],q22[1];
s q22[0];
x q22[1];
h q22[1];
s q22[1];
cx q22[0],q22[1];
sdg q22[0];
h q22[0];
sdg q22[1];
h q22[1];
cx q22[0],q22[1];
z q22[0];
s q22[0];
h q22[0];
s q22[0];
h q22[1];
cx q22[0],q22[1];
sdg q22[0];
h q22[0];
sdg q22[1];
h q22[1];
cx q22[0],q22[1];
h q22[0];
s q22[0];
z q22[0];
sdg q22[0];
h q22[0];
h q22[1];
s q22[1];
x q22[1];
cx q22[0],q22[1];
h q22[0];
s q22[0];
sdg q22[1];
h q22[1];
cx q22[0],q22[1];
z q22[0];
s q22[0];
h q22[1];
s q22[1];
h q22[1];
s q22[1];
cx q22[0],q22[1];
sdg q22[0];
h q22[0];
sdg q22[1];
h q22[1];
cx q22[0],q22[1];
h q22[0];
s q22[0];
h q22[0];
y q22[1];
cx q22[0],q22[1];
h q22[0];
s q22[0];
z q22[0];
sdg q22[1];
h q22[1];
x q22[1];
cx q22[0],q22[1];
sdg q22[0];
h q22[0];
sdg q22[1];
h q22[1];
cx q22[0],q22[1];
x q22[0];
s q22[0];
h q22[1];
s q22[1];
y q22[1];
h q22[1];
cx q22[0],q22[1];
sdg q22[0];
h q22[0];
sdg q22[1];
h q22[1];
cx q22[0],q22[1];
h q22[0];
s q22[0];
h q22[1];
s q22[1];
cx q22[0],q22[1];
y q22[0];
s q22[0];
h q22[0];
s q22[0];
cx q22[0],q22[1];
sdg q22[0];
h q22[0];
sdg q22[1];
h q22[1];
cx q22[0],q22[1];
s q22[0];
x q22[1];
cx q22[0],q22[1];
sdg q22[0];
h q22[0];
x q22[0];
h q22[0];
s q22[0];
h q22[1];
s q22[1];
z q22[1];
cx q22[0],q22[1];
sdg q22[0];
h q22[0];
h q22[0];
cx q22[0],q22[1];
h q22[0];
s q22[0];
sdg q22[1];
h q22[1];
cx q22[0],q22[1];
y q22[0];
h q22[0];
s q22[0];
h q22[1];
s q22[1];
x q22[1];
s q22[1];
h q22[1];
s q22[1];
cx q22[0],q22[1];
x q22[0];
s q22[0];
h q22[0];
s q22[0];
h q22[1];
s q22[1];
cx q22[0],q22[1];
h q22[0];
s q22[0];
sdg q22[1];
h q22[1];
cx q22[0],q22[1];
z q22[0];
h q22[0];
z q22[0];
h q22[0];
z q22[1];
cx q22[0],q22[1];
sdg q22[0];
h q22[0];
h q22[1];
s q22[1];
cx q22[0],q22[1];
s q22[1];
cx q22[0],q22[1];
y q22[0];
sdg q22[1];
h q22[1];
x q22[1];
h q22[1];
cx q22[0],q22[1];
sdg q22[0];
h q22[0];
h q22[1];
s q22[1];
cx q22[0],q22[1];
sdg q22[0];
h q22[0];
y q22[1];
cx q22[0],q22[1];
sdg q22[0];
h q22[0];
sdg q22[1];
h q22[1];
cx q22[0],q22[1];
h q22[0];
s q22[0];
x q22[0];
h q22[1];
s q22[1];
x q22[1];
h q22[1];
s q22[1];
cx q22[0],q22[1];
sdg q22[0];
h q22[0];
s q22[0];
h q22[0];
s q22[0];
sdg q22[1];
h q22[1];
x q22[1];
s q22[1];
h q22[1];
s q22[1];
cx q22[0],q22[1];
sdg q22[0];
h q22[0];
h q22[1];
s q22[1];
cx q22[0],q22[1];
h q22[0];
s q22[0];
x q22[0];
s q22[0];
h q22[0];
s q22[0];
x q22[1];
h q22[1];
s q22[1];
cx q22[0],q22[1];
sdg q22[0];
h q22[0];
h q22[1];
s q22[1];
cx q22[0],q22[1];
h q22[0];
z q22[1];
s q22[1];
cx q22[0],q22[1];
z q22[0];
h q22[0];
s q22[1];
h q22[1];
s q22[1];
cx q22[0],q22[1];
sdg q22[0];
h q22[0];
y q22[0];
h q22[0];
s q22[0];
h q22[1];
s q22[1];
y q22[1];
h q22[1];
cx q22[0],q22[1];
h q22[0];
s q22[0];
sdg q22[1];
h q22[1];
cx q22[0],q22[1];
x q22[0];
z q22[1];
h q22[1];
cx q22[0],q22[1];
sdg q22[0];
h q22[0];
h q22[1];
s q22[1];
cx q22[0],q22[1];
h q22[0];
s q22[0];
y q22[0];
h q22[0];
s q22[0];
x q22[1];
s q22[1];
cx q22[0],q22[1];
sdg q22[0];
h q22[0];
y q22[0];
h q22[0];
x q22[1];
h q22[1];
s q22[1];
cx q22[0],q22[1];
h q22[0];
s q22[0];
z q22[0];
s q22[0];
h q22[0];
s q22[0];
sdg q22[1];
h q22[1];
y q22[1];
s q22[1];
cx q22[0],q22[1];
h q22[0];
s q22[0];
sdg q22[1];
h q22[1];
cx q22[0],q22[1];
z q22[0];
sdg q22[0];
h q22[0];
h q22[1];
s q22[1];
z q22[1];
cx q22[0],q22[1];
sdg q22[0];
h q22[0];
y q22[0];
z q22[1];
cx q22[0],q22[1];
h q22[0];
s q22[0];
sdg q22[1];
h q22[1];
cx q22[0],q22[1];
y q22[0];
s q22[0];
z q22[1];
h q22[1];
s q22[1];
cx q22[0],q22[1];
h q22[0];
s q22[0];
h q22[1];
s q22[1];
cx q22[0],q22[1];
y q22[0];
s q22[0];
h q22[0];
s q22[0];
x q22[1];
h q22[1];
s q22[1];
cx q22[0],q22[1];
sdg q22[0];
h q22[0];
sdg q22[1];
h q22[1];
cx q22[0],q22[1];
z q22[0];
h q22[1];
s q22[1];
z q22[1];
h q22[1];
cx q22[0],q22[1];
sdg q22[0];
h q22[0];
sdg q22[1];
h q22[1];
cx q22[0],q22[1];
s q22[1];
h q22[1];
s q22[1];
cx q22[0],q22[1];
h q22[0];
s q22[0];
sdg q22[0];
h q22[0];
cx q22[0],q22[1];
h q22[0];
s q22[0];
sdg q22[1];
h q22[1];
cx q22[0],q22[1];
z q22[0];
sdg q22[0];
h q22[0];
y q22[1];
s q22[1];
cx q22[0],q22[1];
sdg q22[0];
h q22[0];
y q22[0];
s q22[0];
h q22[1];
s q22[1];
x q22[1];
s q22[1];
cx q22[0],q22[1];
sdg q22[0];
h q22[0];
sdg q22[0];
h q22[0];
sdg q22[1];
h q22[1];
cx q22[0],q22[1];
x q22[0];
h q22[0];
s q22[0];
sdg q22[1];
h q22[1];
sdg q22[1];
h q22[1];
cx q22[0],q22[1];
sdg q22[0];
h q22[0];
sdg q22[1];
h q22[1];
cx q22[0],q22[1];
h q22[0];
s q22[0];
h q22[1];
s q22[1];
cx q22[0],q22[1];
z q22[0];
z q22[1];
cx q22[0],q22[1];
h q22[0];
s q22[0];
sdg q22[1];
h q22[1];
cx q22[0],q22[1];
z q22[0];
h q22[1];
s q22[1];
y q22[1];
h q22[1];
s q22[1];
cx q22[0],q22[1];
sdg q22[0];
h q22[0];
sdg q22[1];
h q22[1];
cx q22[0],q22[1];
h q22[0];
s q22[0];
y q22[1];
s q22[1];
cx q22[0],q22[1];
sdg q22[0];
h q22[0];
h q22[1];
s q22[1];
cx q22[0],q22[1];
h q22[0];
s q22[0];
s q22[0];
x q22[1];
cx q22[0],q22[1];
sdg q22[0];
h q22[0];
sdg q22[1];
h q22[1];
cx q22[0],q22[1];
h q22[0];
s q22[0];
h q22[0];
x q22[1];
s q22[1];
h q22[1];
s q22[1];
cx q22[0],q22[1];
s q22[0];
h q22[0];
s q22[0];
z q22[1];
s q22[1];
h q22[1];
s q22[1];
cx q22[0],q22[1];
sdg q22[0];
h q22[0];
sdg q22[1];
h q22[1];
cx q22[0],q22[1];
h q22[0];
s q22[0];
h q22[1];
s q22[1];
cx q22[0],q22[1];
y q22[0];
s q22[0];
h q22[0];
s q22[0];
y q22[1];
sdg q22[1];
h q22[1];
cx q22[0],q22[1];
h q22[0];
s q22[0];
z q22[0];
sdg q22[1];
h q22[1];
s q22[1];
cx q22[0],q22[1];
sdg q22[0];
h q22[0];
sdg q22[1];
h q22[1];
cx q22[0],q22[1];
h q22[0];
s q22[0];
sdg q22[0];
h q22[0];
y q22[1];
h q22[1];
cx q22[0],q22[1];
sdg q22[0];
h q22[0];
sdg q22[1];
h q22[1];
cx q22[0],q22[1];
h q22[0];
s q22[0];
h q22[1];
s q22[1];
x q22[1];
s q22[1];
cx q22[0],q22[1];
sdg q22[0];
h q22[0];
sdg q22[1];
h q22[1];
cx q22[0],q22[1];
z q22[0];
s q22[0];
h q22[0];
s q22[0];
h q22[1];
s q22[1];
z q22[1];
s q22[1];
cx q22[0],q22[1];
sdg q22[0];
h q22[0];
sdg q22[1];
h q22[1];
cx q22[0],q22[1];
h q22[0];
s q22[0];
h q22[1];
s q22[1];
cx q22[0],q22[1];
x q22[0];
h q22[0];
s q22[0];
y q22[1];
h q22[1];
s q22[1];
cx q22[0],q22[1];
sdg q22[0];
h q22[0];
sdg q22[1];
h q22[1];
cx q22[0],q22[1];
h q22[0];
s q22[0];
y q22[0];
h q22[0];
z q22[1];
s q22[1];
cx q22[0],q22[1];
sdg q22[0];
h q22[0];
h q22[1];
s q22[1];
cx q22[0],q22[1];
h q22[0];
s q22[0];
h q22[0];
x q22[1];
h q22[1];
s q22[1];
cx q22[0],q22[1];
h q22[0];
s q22[0];
sdg q22[1];
h q22[1];
cx q22[0],q22[1];
x q22[0];
y q22[1];
cx q22[0],q22[1];
h q22[0];
s q22[0];
sdg q22[1];
h q22[1];
cx q22[0],q22[1];
x q22[0];
h q22[1];
s q22[1];
x q22[1];
sdg q22[1];
h q22[1];
cx q22[0],q22[1];
sdg q22[0];
h q22[0];
sdg q22[1];
h q22[1];
cx q22[0],q22[1];
h q22[0];
s q22[0];
sdg q22[0];
h q22[0];
h q22[1];
s q22[1];
y q22[1];
sdg q22[1];
h q22[1];
cx q22[0],q22[1];
sdg q22[0];
h q22[0];
sdg q22[1];
h q22[1];
cx q22[0],q22[1];
z q22[0];
h q22[0];
h q22[1];
s q22[1];
y q22[1];
s q22[1];
cx q22[0],q22[1];
sdg q22[0];
h q22[0];
y q22[0];
h q22[1];
s q22[1];
z q22[1];
cx q22[0],q22[1];
sdg q22[0];
h q22[0];
h q22[1];
s q22[1];
cx q22[0],q22[1];
z q22[0];
z q22[1];
h q22[1];
cx q22[0],q22[1];
h q22[0];
s q22[0];
z q22[0];
s q22[0];
h q22[0];
s q22[0];
h q22[1];
s q22[1];
y q22[1];
cx q22[0],q22[1];
sdg q22[0];
h q22[0];
sdg q22[1];
h q22[1];
cx q22[0],q22[1];
h q22[0];
s q22[0];
s q22[0];
h q22[1];
s q22[1];
h q22[1];
s q22[1];
cx q22[0],q22[1];
sdg q22[0];
h q22[0];
sdg q22[1];
h q22[1];
cx q22[0],q22[1];
z q22[0];
s q22[0];
h q22[0];
s q22[0];
s q22[0];
h q22[0];
s q22[0];
h q22[1];
s q22[1];
y q22[1];
s q22[1];
h q22[1];
s q22[1];
cx q22[0],q22[1];
sdg q22[0];
h q22[0];
h q22[1];
s q22[1];
cx q22[0],q22[1];
h q22[0];
s q22[0];
x q22[0];
s q22[0];
s q22[1];
h q22[1];
s q22[1];
cx q22[0],q22[1];
x q22[0];
sdg q22[0];
h q22[0];
z q22[1];
s q22[1];
cx q22[0],q22[1];
sdg q22[0];
h q22[0];
sdg q22[1];
h q22[1];
cx q22[0],q22[1];
h q22[0];
s q22[0];
h q22[1];
s q22[1];
cx q22[0],q22[1];
z q22[0];
h q22[0];
h q22[1];
s q22[1];
cx q22[0],q22[1];
h q22[0];
s q22[0];
y q22[0];
sdg q22[1];
h q22[1];
s q22[1];
cx q22[0],q22[1];
sdg q22[0];
h q22[0];
sdg q22[1];
h q22[1];
cx q22[0],q22[1];
h q22[0];
s q22[0];
z q22[0];
s q22[0];
h q22[0];
s q22[0];
h q22[1];
s q22[1];
z q22[1];
sdg q22[1];
h q22[1];
cx q22[0],q22[1];
sdg q22[0];
h q22[0];
sdg q22[1];
h q22[1];
cx q22[0],q22[1];
z q22[0];
h q22[0];
s q22[0];
y q22[1];
h q22[1];
s q22[1];
cx q22[0],q22[1];
sdg q22[0];
h q22[0];
z q22[0];
s q22[0];
h q22[1];
s q22[1];
sdg q22[1];
h q22[1];
cx q22[0],q22[1];
sdg q22[0];
h q22[0];
s q22[0];
h q22[0];
s q22[0];
h q22[1];
s q22[1];
x q22[1];
s q22[1];
h q22[1];
s q22[1];
cx q22[0],q22[1];
s q22[0];
sdg q22[1];
h q22[1];
cx q22[0],q22[1];
h q22[0];
s q22[0];
sdg q22[1];
h q22[1];
cx q22[0],q22[1];
z q22[0];
h q22[0];
s q22[0];
h q22[1];
s q22[1];
sdg q22[1];
h q22[1];
cx q22[0],q22[1];
h q22[0];
s q22[0];
sdg q22[1];
h q22[1];
cx q22[0],q22[1];
s q22[0];
h q22[0];
s q22[0];
x q22[1];
h q22[1];
cx q22[0],q22[1];
sdg q22[0];
h q22[0];
y q22[0];
s q22[0];
h q22[1];
s q22[1];
z q22[1];
s q22[1];
cx q22[0],q22[1];
sdg q22[0];
h q22[0];
sdg q22[1];
h q22[1];
cx q22[0],q22[1];
h q22[0];
s q22[0];
h q22[1];
s q22[1];
cx q22[0],q22[1];
sdg q22[0];
h q22[0];
y q22[1];
h q22[1];
cx q22[0],q22[1];
h q22[0];
s q22[0];
h q22[0];
z q22[1];
cx q22[0],q22[1];
h q22[0];
s q22[0];
s q22[0];
h q22[0];
s q22[0];
sdg q22[1];
h q22[1];
y q22[1];
sdg q22[1];
h q22[1];
cx q22[0],q22[1];
h q22[0];
s q22[0];
sdg q22[1];
h q22[1];
cx q22[0],q22[1];
y q22[0];
h q22[0];
h q22[1];
s q22[1];
sdg q22[1];
h q22[1];
cx q22[0],q22[1];
sdg q22[0];
h q22[0];
sdg q22[0];
h q22[0];
y q22[1];
sdg q22[1];
h q22[1];
cx q22[0],q22[1];
sdg q22[0];
h q22[0];
sdg q22[1];
h q22[1];
cx q22[0],q22[1];
h q22[0];
s q22[0];
z q22[0];
h q22[0];
s q22[0];
x q22[0];
s q22[0];
h q22[0];
s q22[0];
h q22[1];
s q22[1];
z q22[1];
s q22[1];
x q22[1];
sdg q22[1];
h q22[1];
cx q22[0],q22[1];
x q22[0];
h q22[0];
s q22[0];
sdg q22[1];
h q22[1];
s q22[1];
h q22[1];
s q22[1];
cx q22[0],q22[1];
sdg q22[0];
h q22[0];
h q22[1];
s q22[1];
cx q22[0],q22[1];
y q22[1];
h q22[1];
cx q22[0],q22[1];
sdg q22[0];
h q22[0];
sdg q22[1];
h q22[1];
cx q22[0],q22[1];
z q22[1];
s q22[1];
cx q22[0],q22[1];
sdg q22[0];
h q22[0];
h q22[1];
s q22[1];
cx q22[0],q22[1];
x q22[0];
s q22[0];
h q22[0];
s q22[0];
y q22[1];
h q22[1];
cx q22[0],q22[1];
sdg q22[0];
h q22[0];
h q22[1];
s q22[1];
cx q22[0],q22[1];
x q22[0];
s q22[0];
h q22[0];
s q22[0];
z q22[1];
cx q22[0],q22[1];
h q22[0];
s q22[0];
h q22[1];
s q22[1];
cx q22[0],q22[1];
y q22[0];
s q22[0];
h q22[0];
s q22[0];
z q22[1];
s q22[1];
cx q22[0],q22[1];
h q22[0];
s q22[0];
sdg q22[1];
h q22[1];
cx q22[0],q22[1];
x q22[0];
s q22[0];
h q22[1];
s q22[1];
z q22[1];
s q22[1];
cx q22[0],q22[1];
h q22[0];
s q22[0];
sdg q22[0];
h q22[0];
z q22[0];
s q22[0];
sdg q22[1];
h q22[1];
y q22[1];
s q22[1];
h q22[1];
s q22[1];
x q22[1];
s q22[1];
cx q22[0],q22[1];
sdg q22[0];
h q22[0];
y q22[0];
h q22[0];
s q22[0];
s q22[1];
cx q22[0],q22[1];
sdg q22[0];
h q22[0];
z q22[0];
h q22[1];
s q22[1];
s q22[1];
cx q22[0],q22[1];
sdg q22[0];
h q22[0];
sdg q22[1];
h q22[1];
cx q22[0],q22[1];
h q22[0];
s q22[0];
h q22[0];
h q22[1];
s q22[1];
s q22[1];
cx q22[0],q22[1];
sdg q22[0];
h q22[0];
z q22[0];
sdg q22[1];
h q22[1];
y q22[1];
sdg q22[1];
h q22[1];
cx q22[0],q22[1];
h q22[0];
s q22[0];
s q22[0];
x q22[1];
h q22[1];
s q22[1];
cx q22[0],q22[1];
z q22[0];
h q22[0];
s q22[0];
z q22[1];
sdg q22[1];
h q22[1];
cx q22[0],q22[1];
h q22[0];
s q22[0];
x q22[0];
s q22[0];
h q22[0];
s q22[0];
h q22[1];
s q22[1];
cx q22[0],q22[1];
h q22[0];
s q22[0];
sdg q22[1];
h q22[1];
cx q22[0],q22[1];
x q22[0];
h q22[0];
s q22[0];
h q22[1];
s q22[1];
y q22[1];
h q22[1];
cx q22[0],q22[1];
sdg q22[0];
h q22[0];
sdg q22[1];
h q22[1];
cx q22[0],q22[1];
y q22[0];
h q22[1];
s q22[1];
cx q22[0],q22[1];
h q22[0];
s q22[0];
z q22[0];
h q22[0];
s q22[0];
h q22[0];
s q22[0];
y q22[0];
y q22[1];
z q22[1];
s q22[1];
h q22[1];
s q22[1];
z q22[1];
sdg q22[1];
h q22[1];
cx q22[0],q22[1];
y q22[0];
h q22[0];
h q22[1];
s q22[1];
x q22[1];
h q22[1];
s q22[1];
cx q22[0],q22[1];
h q22[0];
s q22[0];
sdg q22[1];
h q22[1];
cx q22[0],q22[1];
sdg q22[0];
h q22[0];
h q22[1];
s q22[1];
z q22[1];
cx q22[0],q22[1];
x q22[0];
sdg q22[0];
h q22[0];
z q22[1];
cx q22[0],q22[1];
sdg q22[0];
h q22[0];
sdg q22[1];
h q22[1];
cx q22[0],q22[1];
h q22[0];
s q22[0];
h q22[1];
s q22[1];
cx q22[0],q22[1];
z q22[0];
sdg q22[0];
h q22[0];
x q22[1];
s q22[1];
h q22[1];
s q22[1];
cx q22[0],q22[1];
sdg q22[0];
h q22[0];
sdg q22[1];
h q22[1];
cx q22[0],q22[1];
h q22[0];
s q22[0];
h q22[1];
s q22[1];
cx q22[0],q22[1];
y q22[0];
s q22[0];
z q22[1];
h q22[1];
cx q22[0],q22[1];
sdg q22[0];
h q22[0];
s q22[0];
sdg q22[1];
h q22[1];
z q22[1];
h q22[1];
cx q22[0],q22[1];
y q22[0];
sdg q22[0];
h q22[0];
z q22[1];
cx q22[0],q22[1];
sdg q22[0];
h q22[0];
sdg q22[1];
h q22[1];
cx q22[0],q22[1];
z q22[0];
sdg q22[0];
h q22[0];
y q22[0];
s q22[0];
h q22[1];
s q22[1];
z q22[1];
s q22[1];
h q22[1];
s q22[1];
s q22[1];
cx q22[0],q22[1];
sdg q22[0];
h q22[0];
sdg q22[1];
h q22[1];
cx q22[0],q22[1];
h q22[0];
s q22[0];
x q22[0];
s q22[0];
h q22[0];
s q22[0];
x q22[1];
h q22[1];
s q22[1];
cx q22[0],q22[1];
sdg q22[0];
h q22[0];
sdg q22[1];
h q22[1];
cx q22[0],q22[1];
x q22[0];
sdg q22[0];
h q22[0];
h q22[1];
s q22[1];
y q22[1];
s q22[1];
cx q22[0],q22[1];
h q22[0];
s q22[0];
y q22[0];
sdg q22[0];
h q22[0];
sdg q22[1];
h q22[1];
sdg q22[1];
h q22[1];
cx q22[0],q22[1];
x q22[0];
h q22[0];
sdg q22[1];
h q22[1];
z q22[1];
cx q22[0],q22[1];
sdg q22[1];
h q22[1];
x q22[1];
cx q22[0],q22[1];
sdg q22[0];
h q22[0];
y q22[0];
s q22[0];
h q22[0];
s q22[0];
h q22[1];
s q22[1];
s q22[1];
cx q22[0],q22[1];
h q22[0];
s q22[0];
sdg q22[1];
h q22[1];
cx q22[0],q22[1];
y q22[0];
sdg q22[0];
h q22[0];
s q22[1];
cx q22[0],q22[1];
sdg q22[0];
h q22[0];
sdg q22[1];
h q22[1];
cx q22[0],q22[1];
h q22[0];
s q22[0];
h q22[1];
s q22[1];
cx q22[0],q22[1];
y q22[1];
s q22[1];
h q22[1];
s q22[1];
cx q22[0],q22[1];
h q22[0];
s q22[0];
x q22[0];
h q22[1];
s q22[1];
x q22[1];
h q22[1];
s q22[1];
cx q22[0],q22[1];
sdg q22[0];
h q22[0];
z q22[0];
s q22[0];
x q22[1];
s q22[1];
h q22[1];
s q22[1];
cx q22[0],q22[1];
sdg q22[0];
h q22[0];
sdg q22[1];
h q22[1];
cx q22[0],q22[1];
h q22[0];
s q22[0];
h q22[1];
s q22[1];
cx q22[0],q22[1];
x q22[0];
s q22[0];
z q22[1];
s q22[1];
h q22[1];
s q22[1];
cx q22[0],q22[1];
h q22[0];
s q22[0];
sdg q22[1];
h q22[1];
cx q22[0],q22[1];
y q22[0];
sdg q22[0];
h q22[0];
h q22[1];
s q22[1];
x q22[1];
s q22[1];
h q22[1];
s q22[1];
cx q22[0],q22[1];
h q22[0];
s q22[0];
x q22[0];
h q22[0];
s q22[1];
h q22[1];
s q22[1];
cx q22[0],q22[1];
h q22[0];
s q22[0];
s q22[1];
h q22[1];
s q22[1];
cx q22[0],q22[1];
sdg q22[0];
h q22[0];
sdg q22[1];
h q22[1];
cx q22[0],q22[1];
h q22[0];
s q22[0];
h q22[1];
s q22[1];
cx q22[0],q22[1];
y q22[0];
s q22[0];
z q22[1];
s q22[1];
cx q22[0],q22[1];
sdg q22[0];
h q22[0];
sdg q22[1];
h q22[1];
cx q22[0],q22[1];
s q22[0];
h q22[1];
s q22[1];
y q22[1];
cx q22[0],q22[1];
z q22[0];
s q22[0];
sdg q22[1];
h q22[1];
x q22[1];
h q22[1];
cx q22[0],q22[1];
sdg q22[0];
h q22[0];
x q22[0];
s q22[0];
h q22[0];
s q22[0];
sdg q22[1];
h q22[1];
x q22[1];
h q22[1];
cx q22[0],q22[1];
sdg q22[0];
h q22[0];
h q22[1];
s q22[1];
cx q22[0],q22[1];
z q22[0];
s q22[0];
z q22[1];
s q22[1];
h q22[1];
s q22[1];
cx q22[0],q22[1];
sdg q22[0];
h q22[0];
sdg q22[0];
h q22[0];
z q22[1];
h q22[1];
s q22[1];
cx q22[0],q22[1];
sdg q22[0];
h q22[0];
sdg q22[1];
h q22[1];
cx q22[0],q22[1];
h q22[0];
s q22[0];
h q22[1];
s q22[1];
cx q22[0],q22[1];
s q22[0];
h q22[0];
s q22[0];
y q22[1];
h q22[1];
cx q22[0],q22[1];
sdg q22[0];
h q22[0];
sdg q22[1];
h q22[1];
cx q22[0],q22[1];
h q22[0];
s q22[0];
h q22[1];
s q22[1];
cx q22[0],q22[1];
s q22[0];
h q22[0];
s q22[0];
z q22[1];
s q22[1];
cx q22[0],q22[1];
h q22[0];
s q22[0];
sdg q22[1];
h q22[1];
cx q22[0],q22[1];
y q22[0];
h q22[0];
s q22[0];
y q22[1];
s q22[1];
cx q22[0],q22[1];
sdg q22[0];
h q22[0];
sdg q22[1];
h q22[1];
cx q22[0],q22[1];
z q22[0];
h q22[0];
h q22[1];
s q22[1];
z q22[1];
s q22[1];
cx q22[0],q22[1];
h q22[0];
s q22[0];
sdg q22[1];
h q22[1];
cx q22[0],q22[1];
x q22[0];
s q22[0];
h q22[0];
s q22[0];
z q22[1];
s q22[1];
sdg q22[1];
h q22[1];
cx q22[0],q22[1];
sdg q22[0];
h q22[0];
sdg q22[1];
h q22[1];
cx q22[0],q22[1];
y q22[0];
sdg q22[0];
h q22[0];
h q22[1];
s q22[1];
x q22[1];
h q22[1];
s q22[1];
cx q22[0],q22[1];
s q22[0];
sdg q22[1];
h q22[1];
s q22[1];
cx q22[0],q22[1];
sdg q22[0];
h q22[0];
h q22[1];
s q22[1];
cx q22[0],q22[1];
x q22[0];
sdg q22[0];
h q22[0];
x q22[1];
cx q22[0],q22[1];
s q22[1];
cx q22[0],q22[1];
h q22[0];
s q22[0];
z q22[0];
s q22[0];
h q22[1];
s q22[1];
z q22[1];
h q22[1];
cx q22[0],q22[1];
sdg q22[0];
h q22[0];
sdg q22[1];
h q22[1];
cx q22[0],q22[1];
h q22[0];
s q22[0];
h q22[1];
s q22[1];
cx q22[0],q22[1];
s q22[0];
h q22[0];
s q22[0];
y q22[1];
s q22[1];
h q22[1];
s q22[1];
cx q22[0],q22[1];
sdg q22[0];
h q22[0];
sdg q22[1];
h q22[1];
cx q22[0],q22[1];
y q22[0];
h q22[0];
cx q22[0],q22[1];
y q22[0];
h q22[0];
s q22[0];
sdg q22[1];
h q22[1];
h q22[1];
s q22[1];
cx q22[0],q22[1];
h q22[0];
s q22[0];
s q22[0];
h q22[1];
s q22[1];
x q22[1];
h q22[1];
cx q22[0],q22[1];
h q22[0];
s q22[0];
y q22[1];
sdg q22[1];
h q22[1];
cx q22[0],q22[1];
h q22[0];
s q22[0];
sdg q22[1];
h q22[1];
cx q22[0],q22[1];
h q22[0];
x q22[1];
s q22[1];
h q22[1];
s q22[1];
cx q22[0],q22[1];
sdg q22[0];
h q22[0];
sdg q22[1];
h q22[1];
cx q22[0],q22[1];
h q22[0];
s q22[0];
h q22[1];
s q22[1];
cx q22[0],q22[1];
s q22[0];
y q22[1];
s q22[1];
cx q22[0],q22[1];
h q22[0];
s q22[0];
x q22[0];
h q22[0];
s q22[0];
sdg q22[1];
h q22[1];
z q22[1];
cx q22[0],q22[1];
sdg q22[0];
h q22[0];
sdg q22[1];
h q22[1];
cx q22[0],q22[1];
h q22[0];
s q22[0];
h q22[0];
s q22[0];
z q22[1];
s q22[1];
h q22[1];
s q22[1];
cx q22[0],q22[1];
z q22[0];
s q22[0];
x q22[1];
h q22[1];
s q22[1];
cx q22[0],q22[1];
sdg q22[0];
h q22[0];
sdg q22[1];
h q22[1];
