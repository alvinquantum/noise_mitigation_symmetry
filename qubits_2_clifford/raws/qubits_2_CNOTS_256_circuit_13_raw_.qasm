OPENQASM 2.0;
include "qelib1.inc";
qreg q12[2];
y q12[0];
sdg q12[1];
h q12[1];
cx q12[0],q12[1];
h q12[0];
s q12[0];
h q12[1];
s q12[1];
cx q12[0],q12[1];
z q12[0];
h q12[0];
z q12[1];
s q12[1];
h q12[1];
s q12[1];
cx q12[0],q12[1];
h q12[0];
s q12[0];
z q12[0];
h q12[0];
s q12[0];
h q12[1];
s q12[1];
y q12[1];
sdg q12[1];
h q12[1];
cx q12[0],q12[1];
h q12[0];
s q12[0];
h q12[1];
s q12[1];
cx q12[0],q12[1];
z q12[0];
sdg q12[0];
h q12[0];
s q12[1];
cx q12[0],q12[1];
sdg q12[0];
h q12[0];
sdg q12[1];
h q12[1];
cx q12[0],q12[1];
h q12[0];
s q12[1];
h q12[1];
s q12[1];
cx q12[0],q12[1];
sdg q12[0];
h q12[0];
sdg q12[1];
h q12[1];
cx q12[0],q12[1];
h q12[0];
s q12[0];
x q12[0];
h q12[0];
s q12[0];
h q12[1];
s q12[1];
h q12[1];
cx q12[0],q12[1];
sdg q12[0];
h q12[0];
x q12[0];
sdg q12[0];
h q12[0];
sdg q12[1];
h q12[1];
y q12[1];
h q12[1];
s q12[1];
cx q12[0],q12[1];
z q12[0];
h q12[0];
s q12[0];
sdg q12[1];
h q12[1];
x q12[1];
s q12[1];
cx q12[0],q12[1];
h q12[0];
s q12[0];
y q12[0];
sdg q12[1];
h q12[1];
sdg q12[1];
h q12[1];
cx q12[0],q12[1];
h q12[0];
s q12[0];
sdg q12[1];
h q12[1];
cx q12[0],q12[1];
h q12[0];
s q12[0];
z q12[1];
h q12[1];
cx q12[0],q12[1];
sdg q12[0];
h q12[0];
x q12[0];
z q12[1];
sdg q12[1];
h q12[1];
cx q12[0],q12[1];
h q12[0];
s q12[0];
sdg q12[1];
h q12[1];
cx q12[0],q12[1];
z q12[0];
s q12[0];
x q12[1];
s q12[1];
h q12[1];
s q12[1];
cx q12[0],q12[1];
sdg q12[0];
h q12[0];
sdg q12[1];
h q12[1];
cx q12[0],q12[1];
y q12[0];
s q12[0];
h q12[1];
s q12[1];
y q12[1];
s q12[1];
h q12[1];
s q12[1];
cx q12[0],q12[1];
h q12[0];
s q12[0];
x q12[0];
h q12[0];
y q12[1];
sdg q12[1];
h q12[1];
cx q12[0],q12[1];
sdg q12[0];
h q12[0];
sdg q12[1];
h q12[1];
cx q12[0],q12[1];
h q12[0];
s q12[0];
h q12[1];
s q12[1];
cx q12[0],q12[1];
y q12[0];
s q12[0];
h q12[0];
s q12[0];
cx q12[0],q12[1];
sdg q12[0];
h q12[0];
sdg q12[1];
h q12[1];
cx q12[0],q12[1];
x q12[0];
sdg q12[0];
h q12[0];
h q12[1];
s q12[1];
y q12[1];
h q12[1];
s q12[1];
cx q12[0],q12[1];
sdg q12[0];
h q12[0];
sdg q12[1];
h q12[1];
cx q12[0],q12[1];
z q12[0];
h q12[0];
x q12[1];
sdg q12[1];
h q12[1];
cx q12[0],q12[1];
sdg q12[0];
h q12[0];
sdg q12[1];
h q12[1];
cx q12[0],q12[1];
h q12[0];
s q12[0];
h q12[1];
s q12[1];
x q12[1];
h q12[1];
cx q12[0],q12[1];
sdg q12[0];
h q12[0];
sdg q12[1];
h q12[1];
cx q12[0],q12[1];
h q12[0];
s q12[0];
y q12[0];
s q12[0];
h q12[1];
s q12[1];
x q12[1];
sdg q12[1];
h q12[1];
cx q12[0],q12[1];
h q12[0];
s q12[0];
h q12[1];
s q12[1];
cx q12[0],q12[1];
x q12[0];
s q12[0];
y q12[1];
h q12[1];
cx q12[0],q12[1];
h q12[0];
s q12[0];
sdg q12[1];
h q12[1];
cx q12[0],q12[1];
y q12[0];
s q12[0];
h q12[1];
s q12[1];
cx q12[0],q12[1];
h q12[0];
s q12[0];
sdg q12[1];
h q12[1];
cx q12[0],q12[1];
s q12[0];
sdg q12[1];
h q12[1];
cx q12[0],q12[1];
sdg q12[0];
h q12[0];
h q12[1];
s q12[1];
cx q12[0],q12[1];
h q12[0];
s q12[0];
x q12[0];
h q12[1];
cx q12[0],q12[1];
sdg q12[0];
h q12[0];
z q12[0];
h q12[0];
s q12[0];
sdg q12[1];
h q12[1];
y q12[1];
cx q12[0],q12[1];
sdg q12[0];
h q12[0];
sdg q12[1];
h q12[1];
cx q12[0],q12[1];
h q12[0];
x q12[1];
h q12[1];
cx q12[0],q12[1];
sdg q12[0];
h q12[0];
sdg q12[1];
h q12[1];
cx q12[0],q12[1];
h q12[0];
s q12[0];
h q12[1];
s q12[1];
cx q12[0],q12[1];
y q12[0];
sdg q12[0];
h q12[0];
s q12[1];
cx q12[0],q12[1];
sdg q12[0];
h q12[0];
sdg q12[1];
h q12[1];
cx q12[0],q12[1];
h q12[0];
s q12[0];
h q12[1];
s q12[1];
cx q12[0],q12[1];
sdg q12[0];
h q12[0];
z q12[1];
sdg q12[1];
h q12[1];
cx q12[0],q12[1];
h q12[0];
s q12[0];
z q12[0];
s q12[0];
h q12[0];
s q12[0];
y q12[1];
h q12[1];
s q12[1];
cx q12[0],q12[1];
h q12[0];
s q12[0];
sdg q12[1];
h q12[1];
cx q12[0],q12[1];
y q12[0];
y q12[1];
s q12[1];
h q12[1];
s q12[1];
cx q12[0],q12[1];
sdg q12[0];
h q12[0];
x q12[0];
y q12[1];
cx q12[0],q12[1];
h q12[0];
s q12[0];
h q12[1];
s q12[1];
y q12[1];
sdg q12[1];
h q12[1];
cx q12[0],q12[1];
h q12[0];
s q12[0];
z q12[0];
h q12[0];
s q12[0];
x q12[1];
cx q12[0],q12[1];
sdg q12[0];
h q12[0];
h q12[0];
s q12[0];
h q12[1];
s q12[1];
h q12[1];
s q12[1];
cx q12[0],q12[1];
h q12[0];
s q12[0];
x q12[0];
h q12[0];
s q12[0];
z q12[1];
s q12[1];
h q12[1];
s q12[1];
cx q12[0],q12[1];
h q12[0];
s q12[0];
sdg q12[1];
h q12[1];
cx q12[0],q12[1];
h q12[0];
h q12[1];
s q12[1];
x q12[1];
sdg q12[1];
h q12[1];
cx q12[0],q12[1];
z q12[0];
s q12[0];
h q12[0];
s q12[0];
h q12[1];
s q12[1];
y q12[1];
cx q12[0],q12[1];
sdg q12[0];
h q12[0];
h q12[0];
s q12[0];
sdg q12[1];
h q12[1];
s q12[1];
h q12[1];
s q12[1];
cx q12[0],q12[1];
h q12[0];
s q12[0];
z q12[0];
h q12[0];
sdg q12[1];
h q12[1];
z q12[1];
h q12[1];
cx q12[0],q12[1];
sdg q12[0];
h q12[0];
sdg q12[1];
h q12[1];
cx q12[0],q12[1];
h q12[0];
s q12[0];
z q12[0];
sdg q12[0];
h q12[0];
h q12[1];
s q12[1];
z q12[1];
sdg q12[1];
h q12[1];
cx q12[0],q12[1];
sdg q12[0];
h q12[0];
h q12[1];
s q12[1];
cx q12[0],q12[1];
h q12[0];
s q12[0];
h q12[0];
z q12[1];
s q12[1];
cx q12[0],q12[1];
sdg q12[0];
h q12[0];
x q12[0];
sdg q12[0];
h q12[0];
h q12[1];
s q12[1];
x q12[1];
s q12[1];
cx q12[0],q12[1];
h q12[0];
s q12[0];
x q12[0];
sdg q12[0];
h q12[0];
h q12[1];
s q12[1];
y q12[1];
cx q12[0],q12[1];
h q12[0];
s q12[0];
h q12[0];
h q12[1];
s q12[1];
z q12[1];
s q12[1];
cx q12[0],q12[1];
sdg q12[0];
h q12[0];
h q12[0];
s q12[0];
x q12[1];
s q12[1];
cx q12[0],q12[1];
y q12[0];
s q12[0];
h q12[0];
s q12[0];
sdg q12[1];
h q12[1];
y q12[1];
sdg q12[1];
h q12[1];
cx q12[0],q12[1];
h q12[0];
s q12[0];
s q12[0];
h q12[0];
s q12[0];
h q12[1];
s q12[1];
y q12[1];
h q12[1];
cx q12[0],q12[1];
y q12[1];
h q12[1];
s q12[1];
cx q12[0],q12[1];
h q12[0];
s q12[0];
sdg q12[1];
h q12[1];
x q12[1];
sdg q12[1];
h q12[1];
cx q12[0],q12[1];
sdg q12[0];
h q12[0];
h q12[1];
s q12[1];
cx q12[0],q12[1];
z q12[0];
h q12[0];
s q12[0];
y q12[1];
h q12[1];
s q12[1];
cx q12[0],q12[1];
sdg q12[0];
h q12[0];
h q12[0];
s q12[0];
y q12[0];
h q12[0];
sdg q12[1];
h q12[1];
x q12[1];
h q12[1];
s q12[1];
x q12[1];
s q12[1];
h q12[1];
s q12[1];
cx q12[0],q12[1];
h q12[0];
s q12[0];
h q12[1];
s q12[1];
cx q12[0],q12[1];
z q12[0];
s q12[0];
x q12[1];
sdg q12[1];
h q12[1];
z q12[1];
cx q12[0],q12[1];
z q12[0];
h q12[0];
x q12[1];
cx q12[0],q12[1];
sdg q12[0];
h q12[0];
h q12[1];
s q12[1];
cx q12[0],q12[1];
h q12[0];
s q12[0];
s q12[0];
z q12[0];
h q12[0];
sdg q12[1];
h q12[1];
sdg q12[1];
h q12[1];
cx q12[0],q12[1];
h q12[0];
s q12[0];
z q12[0];
h q12[0];
h q12[1];
s q12[1];
x q12[1];
h q12[1];
s q12[1];
cx q12[0],q12[1];
sdg q12[0];
h q12[0];
x q12[0];
h q12[0];
s q12[0];
sdg q12[1];
h q12[1];
z q12[1];
s q12[1];
h q12[1];
s q12[1];
cx q12[0],q12[1];
h q12[0];
s q12[0];
z q12[0];
s q12[0];
sdg q12[1];
h q12[1];
sdg q12[1];
h q12[1];
cx q12[0],q12[1];
sdg q12[0];
h q12[0];
x q12[0];
s q12[0];
cx q12[0],q12[1];
sdg q12[0];
h q12[0];
sdg q12[1];
h q12[1];
cx q12[0],q12[1];
z q12[0];
h q12[0];
h q12[1];
s q12[1];
z q12[1];
cx q12[0],q12[1];
h q12[0];
s q12[0];
z q12[0];
h q12[0];
sdg q12[1];
h q12[1];
h q12[1];
s q12[1];
cx q12[0],q12[1];
h q12[0];
s q12[0];
sdg q12[1];
h q12[1];
cx q12[0],q12[1];
s q12[0];
h q12[0];
s q12[0];
h q12[1];
cx q12[0],q12[1];
sdg q12[0];
h q12[0];
h q12[1];
s q12[1];
x q12[1];
h q12[1];
cx q12[0],q12[1];
sdg q12[0];
h q12[0];
h q12[1];
s q12[1];
cx q12[0],q12[1];
h q12[0];
s q12[0];
x q12[0];
h q12[0];
s q12[1];
cx q12[0],q12[1];
h q12[0];
s q12[0];
sdg q12[1];
h q12[1];
cx q12[0],q12[1];
sdg q12[0];
h q12[0];
x q12[0];
h q12[0];
sdg q12[1];
h q12[1];
cx q12[0],q12[1];
sdg q12[0];
h q12[0];
h q12[1];
s q12[1];
cx q12[0],q12[1];
sdg q12[0];
h q12[0];
y q12[1];
s q12[1];
cx q12[0],q12[1];
sdg q12[0];
h q12[0];
x q12[0];
h q12[0];
sdg q12[1];
h q12[1];
y q12[1];
h q12[1];
s q12[1];
cx q12[0],q12[1];
sdg q12[0];
h q12[0];
sdg q12[1];
h q12[1];
cx q12[0],q12[1];
h q12[0];
s q12[0];
y q12[0];
h q12[0];
h q12[1];
s q12[1];
cx q12[0],q12[1];
x q12[0];
s q12[0];
s q12[1];
h q12[1];
s q12[1];
cx q12[0],q12[1];
sdg q12[0];
h q12[0];
z q12[0];
h q12[0];
s q12[0];
sdg q12[0];
h q12[0];
h q12[1];
s q12[1];
y q12[1];
s q12[1];
z q12[1];
s q12[1];
cx q12[0],q12[1];
h q12[0];
s q12[0];
sdg q12[1];
h q12[1];
cx q12[0],q12[1];
z q12[0];
h q12[0];
s q12[0];
z q12[1];
h q12[1];
cx q12[0],q12[1];
sdg q12[0];
h q12[0];
h q12[0];
sdg q12[1];
h q12[1];
h q12[1];
s q12[1];
cx q12[0],q12[1];
sdg q12[0];
h q12[0];
sdg q12[1];
h q12[1];
cx q12[0],q12[1];
s q12[0];
h q12[0];
s q12[0];
z q12[1];
h q12[1];
cx q12[0],q12[1];
sdg q12[0];
h q12[0];
sdg q12[1];
h q12[1];
cx q12[0],q12[1];
h q12[0];
s q12[0];
y q12[0];
h q12[0];
h q12[1];
s q12[1];
z q12[1];
h q12[1];
s q12[1];
cx q12[0],q12[1];
sdg q12[0];
h q12[0];
x q12[0];
h q12[1];
s q12[1];
x q12[1];
sdg q12[1];
h q12[1];
cx q12[0],q12[1];
sdg q12[0];
h q12[0];
x q12[0];
h q12[1];
s q12[1];
s q12[1];
cx q12[0],q12[1];
sdg q12[0];
h q12[0];
x q12[0];
s q12[0];
sdg q12[1];
h q12[1];
z q12[1];
s q12[1];
cx q12[0],q12[1];
sdg q12[0];
h q12[0];
sdg q12[1];
h q12[1];
cx q12[0],q12[1];
h q12[0];
s q12[0];
h q12[1];
s q12[1];
cx q12[0],q12[1];
y q12[0];
sdg q12[0];
h q12[0];
y q12[1];
h q12[1];
s q12[1];
cx q12[0],q12[1];
sdg q12[0];
h q12[0];
y q12[0];
s q12[0];
sdg q12[1];
h q12[1];
z q12[1];
s q12[1];
h q12[1];
s q12[1];
cx q12[0],q12[1];
sdg q12[0];
h q12[0];
sdg q12[1];
h q12[1];
cx q12[0],q12[1];
h q12[0];
s q12[0];
z q12[0];
s q12[0];
h q12[0];
s q12[0];
y q12[1];
h q12[1];
cx q12[0],q12[1];
sdg q12[0];
h q12[0];
sdg q12[1];
h q12[1];
cx q12[0],q12[1];
h q12[0];
s q12[0];
z q12[0];
s q12[0];
h q12[1];
s q12[1];
h q12[1];
s q12[1];
cx q12[0],q12[1];
h q12[0];
s q12[0];
sdg q12[1];
h q12[1];
z q12[1];
h q12[1];
cx q12[0],q12[1];
sdg q12[0];
h q12[0];
sdg q12[1];
h q12[1];
cx q12[0],q12[1];
h q12[0];
s q12[0];
h q12[1];
s q12[1];
cx q12[0],q12[1];
z q12[0];
s q12[1];
h q12[1];
s q12[1];
cx q12[0],q12[1];
sdg q12[0];
h q12[0];
z q12[0];
s q12[0];
h q12[0];
s q12[0];
sdg q12[1];
h q12[1];
h q12[1];
cx q12[0],q12[1];
h q12[0];
s q12[0];
sdg q12[1];
h q12[1];
x q12[1];
s q12[1];
cx q12[0],q12[1];
x q12[0];
h q12[0];
s q12[0];
h q12[1];
s q12[1];
z q12[1];
sdg q12[1];
h q12[1];
cx q12[0],q12[1];
x q12[0];
s q12[0];
cx q12[0],q12[1];
h q12[0];
s q12[0];
sdg q12[1];
h q12[1];
cx q12[0],q12[1];
z q12[0];
s q12[0];
h q12[1];
s q12[1];
y q12[1];
h q12[1];
s q12[1];
cx q12[0],q12[1];
sdg q12[0];
h q12[0];
sdg q12[1];
h q12[1];
cx q12[0],q12[1];
h q12[0];
s q12[0];
h q12[1];
s q12[1];
cx q12[0],q12[1];
x q12[0];
h q12[0];
z q12[1];
h q12[1];
s q12[1];
cx q12[0],q12[1];
s q12[0];
h q12[0];
s q12[0];
x q12[1];
sdg q12[1];
h q12[1];
cx q12[0],q12[1];
y q12[0];
s q12[0];
sdg q12[1];
h q12[1];
sdg q12[1];
h q12[1];
cx q12[0],q12[1];
sdg q12[0];
h q12[0];
z q12[0];
z q12[0];
s q12[0];
sdg q12[1];
h q12[1];
y q12[1];
h q12[1];
s q12[1];
y q12[1];
s q12[1];
h q12[1];
s q12[1];
cx q12[0],q12[1];
sdg q12[0];
h q12[0];
sdg q12[1];
h q12[1];
cx q12[0],q12[1];
h q12[0];
s q12[0];
h q12[1];
s q12[1];
cx q12[0],q12[1];
x q12[0];
s q12[0];
h q12[0];
s q12[0];
z q12[1];
s q12[1];
h q12[1];
s q12[1];
cx q12[0],q12[1];
h q12[0];
s q12[0];
sdg q12[1];
h q12[1];
cx q12[0],q12[1];
y q12[0];
z q12[1];
cx q12[0],q12[1];
h q12[0];
s q12[0];
z q12[0];
h q12[0];
h q12[1];
cx q12[0],q12[1];
sdg q12[0];
h q12[0];
sdg q12[1];
h q12[1];
cx q12[0],q12[1];
h q12[0];
s q12[0];
y q12[0];
h q12[0];
s q12[0];
h q12[1];
s q12[1];
z q12[1];
cx q12[0],q12[1];
sdg q12[0];
h q12[0];
sdg q12[1];
h q12[1];
cx q12[0],q12[1];
h q12[0];
s q12[0];
z q12[0];
h q12[0];
z q12[1];
sdg q12[1];
h q12[1];
cx q12[0],q12[1];
sdg q12[0];
h q12[0];
z q12[0];
z q12[0];
h q12[0];
s q12[0];
h q12[1];
s q12[1];
y q12[1];
h q12[1];
s q12[1];
cx q12[0],q12[1];
h q12[0];
s q12[0];
h q12[1];
s q12[1];
cx q12[0],q12[1];
x q12[0];
h q12[0];
s q12[0];
x q12[1];
s q12[1];
cx q12[0],q12[1];
h q12[0];
s q12[0];
sdg q12[1];
h q12[1];
cx q12[0],q12[1];
h q12[0];
h q12[1];
s q12[1];
z q12[1];
s q12[1];
cx q12[0],q12[1];
sdg q12[0];
h q12[0];
x q12[0];
sdg q12[0];
h q12[0];
sdg q12[1];
h q12[1];
z q12[1];
sdg q12[1];
h q12[1];
cx q12[0],q12[1];
sdg q12[0];
h q12[0];
sdg q12[1];
h q12[1];
cx q12[0],q12[1];
h q12[0];
s q12[0];
z q12[0];
s q12[0];
h q12[1];
s q12[1];
x q12[1];
sdg q12[1];
h q12[1];
cx q12[0],q12[1];
sdg q12[0];
h q12[0];
y q12[0];
s q12[0];
h q12[0];
s q12[0];
sdg q12[1];
h q12[1];
cx q12[0],q12[1];
sdg q12[0];
h q12[0];
sdg q12[1];
h q12[1];
cx q12[0],q12[1];
h q12[0];
s q12[0];
sdg q12[0];
h q12[0];
h q12[1];
s q12[1];
sdg q12[1];
h q12[1];
cx q12[0],q12[1];
sdg q12[0];
h q12[0];
sdg q12[1];
h q12[1];
cx q12[0],q12[1];
h q12[0];
s q12[0];
h q12[1];
s q12[1];
cx q12[0],q12[1];
z q12[0];
h q12[0];
h q12[0];
x q12[1];
s q12[1];
h q12[1];
s q12[1];
cx q12[0],q12[1];
sdg q12[0];
h q12[0];
h q12[1];
s q12[1];
cx q12[0],q12[1];
s q12[0];
h q12[0];
s q12[0];
z q12[1];
sdg q12[1];
h q12[1];
cx q12[0],q12[1];
z q12[0];
h q12[0];
s q12[0];
s q12[1];
h q12[1];
s q12[1];
cx q12[0],q12[1];
h q12[0];
s q12[0];
sdg q12[1];
h q12[1];
cx q12[0],q12[1];
x q12[0];
h q12[0];
z q12[1];
cx q12[0],q12[1];
sdg q12[0];
h q12[0];
sdg q12[1];
h q12[1];
cx q12[0],q12[1];
h q12[0];
x q12[1];
cx q12[0],q12[1];
x q12[0];
s q12[0];
sdg q12[1];
h q12[1];
s q12[1];
h q12[1];
s q12[1];
cx q12[0],q12[1];
h q12[0];
s q12[0];
sdg q12[1];
h q12[1];
cx q12[0],q12[1];
x q12[0];
h q12[0];
s q12[0];
s q12[1];
h q12[1];
s q12[1];
cx q12[0],q12[1];
sdg q12[0];
h q12[0];
h q12[1];
s q12[1];
cx q12[0],q12[1];
x q12[0];
h q12[0];
x q12[1];
h q12[1];
s q12[1];
cx q12[0],q12[1];
x q12[0];
s q12[0];
h q12[0];
s q12[0];
z q12[0];
z q12[1];
s q12[1];
x q12[1];
s q12[1];
h q12[1];
s q12[1];
cx q12[0],q12[1];
z q12[0];
s q12[0];
h q12[0];
s q12[0];
x q12[1];
cx q12[0],q12[1];
sdg q12[0];
h q12[0];
sdg q12[1];
h q12[1];
cx q12[0],q12[1];
sdg q12[0];
h q12[0];
h q12[1];
s q12[1];
cx q12[0],q12[1];
sdg q12[0];
h q12[0];
h q12[1];
s q12[1];
cx q12[0],q12[1];
s q12[0];
h q12[0];
s q12[0];
s q12[1];
h q12[1];
s q12[1];
cx q12[0],q12[1];
sdg q12[0];
h q12[0];
z q12[0];
s q12[0];
h q12[0];
s q12[0];
sdg q12[1];
h q12[1];
y q12[1];
h q12[1];
s q12[1];
cx q12[0],q12[1];
x q12[0];
h q12[0];
s q12[0];
sdg q12[1];
h q12[1];
y q12[1];
sdg q12[1];
h q12[1];
cx q12[0],q12[1];
x q12[0];
s q12[0];
sdg q12[1];
h q12[1];
y q12[1];
s q12[1];
cx q12[0],q12[1];
h q12[0];
s q12[0];
h q12[1];
s q12[1];
cx q12[0],q12[1];
s q12[0];
h q12[0];
s q12[0];
cx q12[0],q12[1];
sdg q12[0];
h q12[0];
y q12[0];
s q12[0];
h q12[0];
s q12[0];
sdg q12[1];
h q12[1];
x q12[1];
s q12[1];
h q12[1];
s q12[1];
cx q12[0],q12[1];
h q12[0];
s q12[0];
sdg q12[1];
h q12[1];
cx q12[0],q12[1];
y q12[0];
s q12[0];
y q12[1];
s q12[1];
cx q12[0],q12[1];
sdg q12[0];
h q12[0];
sdg q12[1];
h q12[1];
cx q12[0],q12[1];
h q12[0];
s q12[0];
h q12[0];
h q12[1];
s q12[1];
s q12[1];
h q12[1];
s q12[1];
cx q12[0],q12[1];
h q12[0];
s q12[0];
y q12[0];
h q12[0];
s q12[0];
sdg q12[1];
h q12[1];
z q12[1];
s q12[1];
cx q12[0],q12[1];
sdg q12[0];
h q12[0];
sdg q12[1];
h q12[1];
cx q12[0],q12[1];
y q12[0];
h q12[0];
h q12[1];
s q12[1];
y q12[1];
s q12[1];
cx q12[0],q12[1];
h q12[0];
s q12[0];
sdg q12[1];
h q12[1];
cx q12[0],q12[1];
h q12[0];
s q12[0];
h q12[1];
s q12[1];
x q12[1];
h q12[1];
cx q12[0],q12[1];
sdg q12[0];
h q12[0];
h q12[1];
s q12[1];
cx q12[0],q12[1];
y q12[0];
h q12[0];
y q12[1];
h q12[1];
cx q12[0],q12[1];
sdg q12[0];
h q12[0];
x q12[0];
s q12[0];
h q12[0];
s q12[0];
h q12[1];
s q12[1];
s q12[1];
cx q12[0],q12[1];
sdg q12[0];
h q12[0];
x q12[0];
h q12[0];
s q12[0];
sdg q12[1];
h q12[1];
s q12[1];
cx q12[0],q12[1];
sdg q12[0];
h q12[0];
h q12[1];
s q12[1];
cx q12[0],q12[1];
h q12[0];
s q12[0];
h q12[0];
s q12[0];
z q12[1];
s q12[1];
h q12[1];
s q12[1];
cx q12[0],q12[1];
h q12[0];
s q12[0];
x q12[0];
sdg q12[0];
h q12[0];
h q12[1];
s q12[1];
z q12[1];
h q12[1];
cx q12[0],q12[1];
sdg q12[0];
h q12[0];
sdg q12[1];
h q12[1];
cx q12[0],q12[1];
h q12[0];
s q12[0];
y q12[0];
s q12[0];
x q12[0];
s q12[0];
h q12[0];
s q12[0];
h q12[1];
s q12[1];
sdg q12[1];
h q12[1];
cx q12[0],q12[1];
sdg q12[0];
h q12[0];
sdg q12[1];
h q12[1];
cx q12[0],q12[1];
s q12[0];
h q12[0];
s q12[0];
h q12[1];
s q12[1];
sdg q12[1];
h q12[1];
cx q12[0],q12[1];
sdg q12[0];
h q12[0];
h q12[1];
s q12[1];
cx q12[0],q12[1];
x q12[0];
s q12[0];
h q12[0];
s q12[0];
z q12[1];
h q12[1];
s q12[1];
cx q12[0],q12[1];
h q12[0];
s q12[0];
sdg q12[0];
h q12[0];
z q12[1];
s q12[1];
cx q12[0],q12[1];
sdg q12[0];
h q12[0];
x q12[1];
h q12[1];
cx q12[0],q12[1];
h q12[0];
s q12[0];
sdg q12[1];
h q12[1];
cx q12[0],q12[1];
sdg q12[0];
h q12[0];
x q12[1];
h q12[1];
s q12[1];
cx q12[0],q12[1];
sdg q12[0];
h q12[0];
s q12[0];
h q12[0];
s q12[0];
y q12[1];
sdg q12[1];
h q12[1];
cx q12[0],q12[1];
x q12[0];
sdg q12[0];
h q12[0];
h q12[1];
s q12[1];
z q12[1];
s q12[1];
h q12[1];
s q12[1];
cx q12[0],q12[1];
z q12[0];
sdg q12[0];
h q12[0];
h q12[0];
h q12[1];
s q12[1];
z q12[1];
sdg q12[1];
h q12[1];
h q12[1];
cx q12[0],q12[1];
sdg q12[0];
h q12[0];
h q12[1];
s q12[1];
cx q12[0],q12[1];
y q12[0];
sdg q12[0];
h q12[0];
s q12[1];
cx q12[0],q12[1];
sdg q12[0];
h q12[0];
sdg q12[1];
h q12[1];
cx q12[0],q12[1];
h q12[0];
s q12[0];
z q12[0];
h q12[0];
h q12[1];
s q12[1];
x q12[1];
h q12[1];
s q12[1];
cx q12[0],q12[1];
sdg q12[0];
h q12[0];
h q12[1];
s q12[1];
cx q12[0],q12[1];
h q12[0];
s q12[0];
z q12[0];
h q12[0];
s q12[0];
z q12[1];
h q12[1];
s q12[1];
cx q12[0],q12[1];
h q12[0];
s q12[0];
h q12[1];
s q12[1];
cx q12[0],q12[1];
sdg q12[0];
h q12[0];
x q12[1];
sdg q12[1];
h q12[1];
cx q12[0],q12[1];
h q12[0];
s q12[0];
x q12[0];
sdg q12[0];
h q12[0];
h q12[1];
s q12[1];
y q12[1];
s q12[1];
cx q12[0],q12[1];
sdg q12[0];
h q12[0];
sdg q12[1];
h q12[1];
cx q12[0],q12[1];
h q12[0];
s q12[0];
z q12[0];
h q12[0];
x q12[1];
cx q12[0],q12[1];
sdg q12[0];
h q12[0];
x q12[0];
sdg q12[0];
h q12[0];
sdg q12[1];
h q12[1];
z q12[1];
s q12[1];
cx q12[0],q12[1];
sdg q12[0];
h q12[0];
z q12[0];
sdg q12[0];
h q12[0];
y q12[1];
h q12[1];
s q12[1];
cx q12[0],q12[1];
h q12[0];
s q12[0];
sdg q12[1];
h q12[1];
cx q12[0],q12[1];
sdg q12[0];
h q12[0];
h q12[1];
s q12[1];
y q12[1];
cx q12[0],q12[1];
sdg q12[0];
h q12[0];
sdg q12[1];
h q12[1];
cx q12[0],q12[1];
h q12[0];
s q12[0];
y q12[0];
s q12[0];
h q12[0];
s q12[0];
h q12[1];
s q12[1];
x q12[1];
s q12[1];
h q12[1];
s q12[1];
cx q12[0],q12[1];
sdg q12[0];
h q12[0];
sdg q12[1];
h q12[1];
cx q12[0],q12[1];
h q12[0];
s q12[0];
y q12[0];
s q12[0];
h q12[1];
s q12[1];
x q12[1];
s q12[1];
h q12[1];
s q12[1];
cx q12[0],q12[1];
h q12[0];
s q12[0];
sdg q12[1];
h q12[1];
cx q12[0],q12[1];
sdg q12[0];
h q12[0];
z q12[1];
h q12[1];
cx q12[0],q12[1];
h q12[0];
s q12[0];
h q12[1];
s q12[1];
z q12[1];
h q12[1];
s q12[1];
cx q12[0],q12[1];
sdg q12[0];
h q12[0];
sdg q12[1];
h q12[1];
cx q12[0],q12[1];
h q12[0];
s q12[0];
h q12[1];
s q12[1];
cx q12[0],q12[1];
y q12[0];
h q12[0];
s q12[0];
h q12[1];
cx q12[0],q12[1];
sdg q12[0];
h q12[0];
sdg q12[1];
h q12[1];
cx q12[0],q12[1];
h q12[0];
s q12[0];
y q12[0];
h q12[0];
s q12[0];
y q12[1];
cx q12[0],q12[1];
h q12[0];
s q12[0];
sdg q12[1];
h q12[1];
cx q12[0],q12[1];
x q12[0];
h q12[1];
s q12[1];
cx q12[0],q12[1];
h q12[0];
s q12[0];
h q12[1];
s q12[1];
cx q12[0],q12[1];
z q12[0];
h q12[0];
s q12[0];
x q12[1];
h q12[1];
s q12[1];
cx q12[0],q12[1];
h q12[0];
s q12[0];
sdg q12[1];
h q12[1];
cx q12[0],q12[1];
y q12[0];
s q12[0];
h q12[1];
s q12[1];
sdg q12[1];
h q12[1];
cx q12[0],q12[1];
sdg q12[0];
h q12[0];
sdg q12[1];
h q12[1];
cx q12[0],q12[1];
h q12[0];
s q12[0];
h q12[0];
y q12[1];
h q12[1];
cx q12[0],q12[1];
sdg q12[0];
h q12[0];
sdg q12[1];
h q12[1];
cx q12[0],q12[1];
y q12[0];
s q12[0];
z q12[1];
sdg q12[1];
h q12[1];
cx q12[0],q12[1];
h q12[0];
s q12[0];
sdg q12[1];
h q12[1];
cx q12[0],q12[1];
s q12[0];
y q12[1];
sdg q12[1];
h q12[1];
cx q12[0],q12[1];
h q12[0];
s q12[0];
h q12[1];
s q12[1];
cx q12[0],q12[1];
sdg q12[0];
h q12[0];
sdg q12[1];
h q12[1];
cx q12[0],q12[1];
x q12[0];
s q12[0];
z q12[1];
h q12[1];
s q12[1];
cx q12[0],q12[1];
sdg q12[0];
h q12[0];
sdg q12[1];
h q12[1];
