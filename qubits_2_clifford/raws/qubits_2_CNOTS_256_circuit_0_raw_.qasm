OPENQASM 2.0;
include "qelib1.inc";
qreg q0[2];
y q0[0];
s q0[0];
h q0[0];
s q0[0];
x q0[1];
s q0[1];
h q0[1];
s q0[1];
cx q0[0],q0[1];
sdg q0[0];
h q0[0];
y q0[0];
z q0[1];
s q0[1];
cx q0[0],q0[1];
x q0[0];
sdg q0[1];
h q0[1];
z q0[1];
h q0[1];
s q0[1];
cx q0[0],q0[1];
sdg q0[0];
h q0[0];
sdg q0[1];
h q0[1];
cx q0[0],q0[1];
x q0[0];
sdg q0[0];
h q0[0];
h q0[1];
s q0[1];
y q0[1];
s q0[1];
cx q0[0],q0[1];
sdg q0[0];
h q0[0];
sdg q0[1];
h q0[1];
cx q0[0],q0[1];
h q0[0];
s q0[0];
h q0[1];
s q0[1];
cx q0[0],q0[1];
z q0[0];
y q0[1];
h q0[1];
cx q0[0],q0[1];
sdg q0[0];
h q0[0];
h q0[1];
s q0[1];
cx q0[0],q0[1];
h q0[0];
s q0[0];
cx q0[0],q0[1];
h q0[0];
s q0[0];
z q0[0];
s q0[0];
h q0[0];
s q0[0];
y q0[0];
h q0[0];
s q0[0];
sdg q0[1];
h q0[1];
x q0[1];
z q0[1];
s q0[1];
h q0[1];
s q0[1];
cx q0[0],q0[1];
sdg q0[0];
h q0[0];
h q0[1];
s q0[1];
cx q0[0],q0[1];
h q0[0];
s q0[0];
y q0[0];
h q0[0];
s q0[0];
x q0[1];
s q0[1];
h q0[1];
s q0[1];
cx q0[0],q0[1];
sdg q0[0];
h q0[0];
h q0[1];
s q0[1];
cx q0[0],q0[1];
h q0[0];
s q0[0];
z q0[0];
x q0[1];
cx q0[0],q0[1];
sdg q0[0];
h q0[0];
sdg q0[1];
h q0[1];
cx q0[0],q0[1];
h q0[0];
s q0[0];
y q0[1];
sdg q0[1];
h q0[1];
cx q0[0],q0[1];
h q0[0];
s q0[0];
sdg q0[1];
h q0[1];
cx q0[0],q0[1];
y q0[0];
s q0[0];
x q0[1];
h q0[1];
cx q0[0],q0[1];
sdg q0[0];
h q0[0];
sdg q0[1];
h q0[1];
cx q0[0],q0[1];
z q0[0];
h q0[0];
s q0[0];
h q0[1];
s q0[1];
x q0[1];
sdg q0[1];
h q0[1];
cx q0[0],q0[1];
h q0[0];
s q0[0];
y q0[0];
s q0[0];
h q0[0];
s q0[0];
sdg q0[1];
h q0[1];
s q0[1];
cx q0[0],q0[1];
sdg q0[0];
h q0[0];
h q0[1];
s q0[1];
cx q0[0],q0[1];
z q0[0];
h q0[0];
s q0[0];
h q0[1];
cx q0[0],q0[1];
h q0[0];
s q0[0];
h q0[1];
s q0[1];
cx q0[0],q0[1];
x q0[0];
s q0[1];
h q0[1];
s q0[1];
cx q0[0],q0[1];
sdg q0[0];
h q0[0];
sdg q0[1];
h q0[1];
cx q0[0],q0[1];
y q0[0];
sdg q0[0];
h q0[0];
h q0[1];
s q0[1];
x q0[1];
s q0[1];
h q0[1];
s q0[1];
cx q0[0],q0[1];
sdg q0[0];
h q0[0];
h q0[1];
s q0[1];
y q0[1];
h q0[1];
s q0[1];
cx q0[0],q0[1];
h q0[0];
s q0[0];
h q0[0];
sdg q0[1];
h q0[1];
z q0[1];
s q0[1];
cx q0[0],q0[1];
h q0[0];
s q0[0];
z q0[0];
s q0[0];
z q0[1];
s q0[1];
h q0[1];
s q0[1];
cx q0[0],q0[1];
h q0[0];
s q0[0];
sdg q0[1];
h q0[1];
cx q0[0],q0[1];
y q0[0];
h q0[0];
s q0[0];
h q0[1];
s q0[1];
x q0[1];
sdg q0[1];
h q0[1];
cx q0[0],q0[1];
sdg q0[0];
h q0[0];
y q0[0];
s q0[0];
z q0[1];
h q0[1];
cx q0[0],q0[1];
sdg q0[0];
h q0[0];
z q0[0];
h q0[0];
s q0[0];
h q0[1];
s q0[1];
y q0[1];
cx q0[0],q0[1];
h q0[0];
s q0[0];
x q0[0];
s q0[0];
s q0[1];
h q0[1];
s q0[1];
cx q0[0],q0[1];
h q0[0];
s q0[0];
z q0[0];
s q0[0];
z q0[1];
s q0[1];
cx q0[0],q0[1];
sdg q0[0];
h q0[0];
sdg q0[1];
h q0[1];
cx q0[0],q0[1];
z q0[0];
y q0[1];
h q0[1];
cx q0[0],q0[1];
sdg q0[0];
h q0[0];
sdg q0[1];
h q0[1];
cx q0[0],q0[1];
h q0[0];
s q0[0];
s q0[0];
h q0[0];
s q0[0];
h q0[1];
s q0[1];
x q0[1];
h q0[1];
cx q0[0],q0[1];
sdg q0[0];
h q0[0];
h q0[1];
s q0[1];
cx q0[0],q0[1];
y q0[0];
sdg q0[0];
h q0[0];
y q0[1];
sdg q0[1];
h q0[1];
cx q0[0],q0[1];
sdg q0[0];
h q0[0];
sdg q0[1];
h q0[1];
cx q0[0],q0[1];
y q0[0];
s q0[0];
y q0[1];
cx q0[0],q0[1];
sdg q0[0];
h q0[0];
h q0[0];
s q0[0];
sdg q0[1];
h q0[1];
x q0[1];
cx q0[0],q0[1];
h q0[0];
s q0[0];
sdg q0[1];
h q0[1];
cx q0[0],q0[1];
x q0[0];
sdg q0[0];
h q0[0];
y q0[1];
cx q0[0],q0[1];
sdg q0[0];
h q0[0];
sdg q0[1];
h q0[1];
cx q0[0],q0[1];
y q0[0];
sdg q0[0];
h q0[0];
sdg q0[1];
h q0[1];
cx q0[0],q0[1];
h q0[0];
s q0[0];
sdg q0[1];
h q0[1];
cx q0[0],q0[1];
y q0[0];
s q0[0];
x q0[1];
s q0[1];
h q0[1];
s q0[1];
cx q0[0],q0[1];
sdg q0[0];
h q0[0];
z q0[0];
h q0[0];
sdg q0[1];
h q0[1];
x q0[1];
h q0[1];
cx q0[0],q0[1];
h q0[0];
s q0[0];
y q0[0];
s q0[0];
h q0[0];
s q0[0];
z q0[1];
cx q0[0],q0[1];
h q0[0];
s q0[0];
h q0[1];
s q0[1];
cx q0[0],q0[1];
h q0[0];
z q0[1];
cx q0[0],q0[1];
h q0[0];
s q0[0];
y q0[1];
sdg q0[1];
h q0[1];
cx q0[0],q0[1];
y q0[0];
s q0[0];
h q0[0];
s q0[0];
h q0[1];
s q0[1];
x q0[1];
sdg q0[1];
h q0[1];
cx q0[0],q0[1];
sdg q0[0];
h q0[0];
h q0[1];
s q0[1];
cx q0[0],q0[1];
h q0[0];
s q0[0];
h q0[0];
s q0[0];
z q0[1];
h q0[1];
cx q0[0],q0[1];
sdg q0[0];
h q0[0];
x q0[0];
s q0[0];
sdg q0[1];
h q0[1];
x q0[1];
sdg q0[1];
h q0[1];
cx q0[0],q0[1];
x q0[0];
h q0[1];
s q0[1];
x q0[1];
h q0[1];
cx q0[0],q0[1];
sdg q0[0];
h q0[0];
h q0[1];
s q0[1];
cx q0[0],q0[1];
y q0[0];
h q0[0];
s q0[0];
x q0[1];
s q0[1];
cx q0[0],q0[1];
sdg q0[0];
h q0[0];
sdg q0[1];
h q0[1];
cx q0[0],q0[1];
h q0[0];
s q0[0];
h q0[1];
s q0[1];
y q0[1];
s q0[1];
h q0[1];
s q0[1];
cx q0[0],q0[1];
h q0[0];
s q0[0];
sdg q0[1];
h q0[1];
cx q0[0],q0[1];
s q0[0];
x q0[1];
sdg q0[1];
h q0[1];
cx q0[0],q0[1];
h q0[0];
s q0[0];
h q0[0];
sdg q0[1];
h q0[1];
z q0[1];
s q0[1];
cx q0[0],q0[1];
h q0[0];
s q0[0];
h q0[1];
s q0[1];
cx q0[0],q0[1];
x q0[0];
h q0[0];
s q0[0];
x q0[1];
s q0[1];
cx q0[0],q0[1];
sdg q0[0];
h q0[0];
z q0[0];
h q0[0];
sdg q0[1];
h q0[1];
z q0[1];
h q0[1];
cx q0[0],q0[1];
h q0[0];
s q0[0];
x q0[0];
sdg q0[0];
h q0[0];
z q0[0];
sdg q0[1];
h q0[1];
z q0[1];
h q0[1];
s q0[1];
x q0[1];
s q0[1];
cx q0[0],q0[1];
y q0[0];
sdg q0[0];
h q0[0];
h q0[1];
s q0[1];
y q0[1];
cx q0[0],q0[1];
x q0[0];
h q0[0];
h q0[1];
s q0[1];
x q0[1];
sdg q0[1];
h q0[1];
cx q0[0],q0[1];
h q0[0];
s q0[0];
z q0[0];
sdg q0[0];
h q0[0];
x q0[1];
s q0[1];
cx q0[0],q0[1];
sdg q0[0];
h q0[0];
z q0[0];
h q0[1];
s q0[1];
cx q0[0],q0[1];
sdg q0[0];
h q0[0];
z q0[0];
h q0[1];
s q0[1];
s q0[1];
h q0[1];
s q0[1];
cx q0[0],q0[1];
sdg q0[0];
h q0[0];
h q0[1];
s q0[1];
cx q0[0],q0[1];
z q0[1];
s q0[1];
h q0[1];
s q0[1];
cx q0[0],q0[1];
h q0[0];
s q0[0];
s q0[0];
h q0[0];
s q0[0];
h q0[1];
s q0[1];
h q0[1];
cx q0[0],q0[1];
x q0[0];
h q0[0];
s q0[0];
h q0[1];
s q0[1];
sdg q0[1];
h q0[1];
cx q0[0],q0[1];
h q0[0];
s q0[0];
sdg q0[0];
h q0[0];
sdg q0[1];
h q0[1];
s q0[1];
h q0[1];
s q0[1];
cx q0[0],q0[1];
sdg q0[0];
h q0[0];
x q0[0];
x q0[1];
h q0[1];
s q0[1];
cx q0[0],q0[1];
h q0[0];
s q0[0];
z q0[0];
s q0[0];
h q0[0];
s q0[0];
y q0[1];
s q0[1];
h q0[1];
s q0[1];
cx q0[0],q0[1];
h q0[0];
s q0[0];
h q0[1];
s q0[1];
cx q0[0],q0[1];
s q0[0];
h q0[0];
s q0[0];
z q0[1];
sdg q0[1];
h q0[1];
cx q0[0],q0[1];
sdg q0[0];
h q0[0];
x q0[0];
h q0[0];
sdg q0[1];
h q0[1];
z q0[1];
sdg q0[1];
h q0[1];
cx q0[0],q0[1];
h q0[0];
s q0[0];
x q0[1];
cx q0[0],q0[1];
h q0[0];
s q0[0];
sdg q0[0];
h q0[0];
h q0[1];
s q0[1];
y q0[1];
cx q0[0],q0[1];
y q0[0];
s q0[0];
sdg q0[1];
h q0[1];
x q0[1];
s q0[1];
cx q0[0],q0[1];
sdg q0[0];
h q0[0];
h q0[1];
s q0[1];
cx q0[0],q0[1];
x q0[0];
s q0[0];
h q0[0];
s q0[0];
x q0[1];
h q0[1];
s q0[1];
cx q0[0],q0[1];
sdg q0[0];
h q0[0];
z q0[0];
sdg q0[0];
h q0[0];
sdg q0[1];
h q0[1];
x q0[1];
cx q0[0],q0[1];
sdg q0[0];
h q0[0];
sdg q0[1];
h q0[1];
cx q0[0],q0[1];
h q0[0];
s q0[0];
y q0[0];
s q0[0];
h q0[1];
s q0[1];
x q0[1];
h q0[1];
s q0[1];
cx q0[0],q0[1];
sdg q0[0];
h q0[0];
sdg q0[1];
h q0[1];
cx q0[0],q0[1];
h q0[0];
s q0[0];
z q0[0];
h q0[0];
s q0[0];
h q0[1];
s q0[1];
x q0[1];
s q0[1];
cx q0[0],q0[1];
sdg q0[0];
h q0[0];
z q0[0];
s q0[0];
h q0[0];
s q0[0];
h q0[1];
s q0[1];
z q0[1];
cx q0[0],q0[1];
h q0[0];
s q0[0];
h q0[1];
s q0[1];
cx q0[0],q0[1];
z q0[0];
s q0[0];
h q0[0];
s q0[0];
x q0[1];
sdg q0[1];
h q0[1];
cx q0[0],q0[1];
sdg q0[0];
h q0[0];
x q0[0];
s q0[0];
sdg q0[1];
h q0[1];
y q0[1];
h q0[1];
s q0[1];
cx q0[0],q0[1];
sdg q0[0];
h q0[0];
sdg q0[1];
h q0[1];
cx q0[0],q0[1];
h q0[0];
s q0[0];
x q0[0];
s q0[0];
x q0[1];
h q0[1];
s q0[1];
cx q0[0],q0[1];
y q0[0];
h q0[0];
s q0[1];
cx q0[0],q0[1];
h q0[0];
s q0[0];
sdg q0[1];
h q0[1];
cx q0[0],q0[1];
y q0[0];
s q0[0];
h q0[0];
s q0[0];
h q0[1];
s q0[1];
y q0[1];
h q0[1];
cx q0[0],q0[1];
h q0[0];
s q0[0];
h q0[1];
s q0[1];
cx q0[0],q0[1];
x q0[0];
s q0[0];
h q0[0];
s q0[0];
x q0[1];
h q0[1];
s q0[1];
cx q0[0],q0[1];
h q0[0];
s q0[0];
x q0[0];
h q0[0];
sdg q0[1];
h q0[1];
x q0[1];
sdg q0[1];
h q0[1];
cx q0[0],q0[1];
sdg q0[0];
h q0[0];
z q0[0];
s q0[0];
h q0[0];
s q0[0];
x q0[1];
cx q0[0],q0[1];
sdg q0[0];
h q0[0];
sdg q0[1];
h q0[1];
cx q0[0],q0[1];
h q0[0];
s q0[0];
h q0[1];
s q0[1];
cx q0[0],q0[1];
s q0[0];
z q0[1];
h q0[1];
s q0[1];
cx q0[0],q0[1];
sdg q0[0];
h q0[0];
z q0[0];
s q0[0];
x q0[0];
h q0[0];
h q0[1];
s q0[1];
y q0[1];
s q0[1];
cx q0[0],q0[1];
sdg q0[0];
h q0[0];
sdg q0[1];
h q0[1];
cx q0[0],q0[1];
h q0[0];
s q0[0];
y q0[0];
sdg q0[0];
h q0[0];
sdg q0[1];
h q0[1];
cx q0[0],q0[1];
sdg q0[0];
h q0[0];
x q0[0];
sdg q0[1];
h q0[1];
z q0[1];
h q0[1];
cx q0[0],q0[1];
h q0[0];
s q0[0];
y q0[0];
s q0[0];
h q0[0];
s q0[0];
s q0[1];
h q0[1];
s q0[1];
cx q0[0],q0[1];
h q0[0];
s q0[0];
sdg q0[1];
h q0[1];
cx q0[0],q0[1];
z q0[0];
s q0[0];
z q0[1];
h q0[1];
cx q0[0],q0[1];
y q0[0];
sdg q0[0];
h q0[0];
h q0[1];
s q0[1];
x q0[1];
s q0[1];
h q0[1];
s q0[1];
cx q0[0],q0[1];
h q0[0];
s q0[0];
sdg q0[1];
h q0[1];
cx q0[0],q0[1];
z q0[0];
h q0[0];
sdg q0[1];
h q0[1];
cx q0[0],q0[1];
h q0[0];
s q0[0];
h q0[1];
s q0[1];
cx q0[0],q0[1];
y q0[0];
h q0[0];
cx q0[0],q0[1];
s q0[0];
h q0[0];
s q0[0];
y q0[1];
s q0[1];
cx q0[0],q0[1];
sdg q0[0];
h q0[0];
h q0[1];
s q0[1];
cx q0[0],q0[1];
z q0[0];
s q0[0];
h q0[0];
s q0[0];
x q0[1];
h q0[1];
s q0[1];
cx q0[0],q0[1];
h q0[0];
s q0[0];
sdg q0[1];
h q0[1];
cx q0[0],q0[1];
x q0[0];
h q0[0];
h q0[1];
s q0[1];
y q0[1];
s q0[1];
cx q0[0],q0[1];
sdg q0[0];
h q0[0];
h q0[1];
s q0[1];
cx q0[0],q0[1];
h q0[0];
s q0[0];
z q0[0];
sdg q0[0];
h q0[0];
y q0[1];
h q0[1];
s q0[1];
cx q0[0],q0[1];
z q0[0];
s q0[0];
sdg q0[1];
h q0[1];
sdg q0[1];
h q0[1];
cx q0[0],q0[1];
sdg q0[0];
h q0[0];
sdg q0[1];
h q0[1];
cx q0[0],q0[1];
h q0[0];
s q0[0];
h q0[1];
s q0[1];
cx q0[0],q0[1];
s q0[0];
h q0[0];
s q0[0];
z q0[1];
h q0[1];
cx q0[0],q0[1];
sdg q0[0];
h q0[0];
s q0[0];
h q0[0];
s q0[0];
h q0[1];
s q0[1];
y q0[1];
s q0[1];
h q0[1];
s q0[1];
cx q0[0],q0[1];
s q0[0];
sdg q0[1];
h q0[1];
x q0[1];
h q0[1];
s q0[1];
cx q0[0],q0[1];
h q0[0];
s q0[0];
h q0[0];
s q0[0];
z q0[1];
sdg q0[1];
h q0[1];
cx q0[0],q0[1];
sdg q0[0];
h q0[0];
sdg q0[0];
h q0[0];
x q0[1];
sdg q0[1];
h q0[1];
cx q0[0],q0[1];
h q0[0];
s q0[0];
z q0[0];
h q0[0];
s q0[0];
h q0[1];
s q0[1];
z q0[1];
sdg q0[1];
h q0[1];
cx q0[0],q0[1];
sdg q0[0];
h q0[0];
sdg q0[1];
h q0[1];
cx q0[0],q0[1];
z q0[0];
h q0[1];
s q0[1];
y q0[1];
sdg q0[1];
h q0[1];
cx q0[0],q0[1];
h q0[0];
z q0[1];
h q0[1];
cx q0[0],q0[1];
sdg q0[0];
h q0[0];
sdg q0[1];
h q0[1];
cx q0[0],q0[1];
x q0[0];
h q0[0];
s q0[0];
h q0[1];
s q0[1];
z q0[1];
sdg q0[1];
h q0[1];
cx q0[0],q0[1];
sdg q0[0];
h q0[0];
h q0[1];
s q0[1];
cx q0[0],q0[1];
h q0[0];
s q0[0];
y q0[0];
s q0[0];
h q0[0];
s q0[0];
y q0[1];
cx q0[0],q0[1];
h q0[0];
s q0[0];
x q0[0];
z q0[1];
sdg q0[1];
h q0[1];
cx q0[0],q0[1];
sdg q0[0];
h q0[0];
sdg q0[1];
h q0[1];
cx q0[0],q0[1];
x q0[0];
h q0[0];
h q0[1];
s q0[1];
z q0[1];
cx q0[0],q0[1];
h q0[0];
s q0[0];
x q0[0];
sdg q0[0];
h q0[0];
y q0[1];
h q0[1];
cx q0[0],q0[1];
sdg q0[0];
h q0[0];
h q0[1];
s q0[1];
cx q0[0],q0[1];
h q0[0];
s q0[0];
z q0[0];
sdg q0[0];
h q0[0];
sdg q0[1];
h q0[1];
cx q0[0],q0[1];
sdg q0[0];
h q0[0];
sdg q0[1];
h q0[1];
cx q0[0],q0[1];
x q0[0];
h q0[0];
s q0[0];
h q0[1];
s q0[1];
z q0[1];
sdg q0[1];
h q0[1];
cx q0[0],q0[1];
sdg q0[0];
h q0[0];
sdg q0[1];
h q0[1];
cx q0[0],q0[1];
x q0[0];
z q0[1];
h q0[1];
cx q0[0],q0[1];
h q0[0];
s q0[0];
sdg q0[1];
h q0[1];
cx q0[0],q0[1];
h q0[0];
s q0[0];
z q0[1];
h q0[1];
cx q0[0],q0[1];
y q0[0];
h q0[0];
s q0[0];
y q0[1];
cx q0[0],q0[1];
sdg q0[0];
h q0[0];
h q0[1];
s q0[1];
cx q0[0],q0[1];
h q0[0];
s q0[0];
s q0[0];
h q0[1];
s q0[1];
cx q0[0],q0[1];
sdg q0[0];
h q0[0];
sdg q0[1];
h q0[1];
cx q0[0],q0[1];
h q0[0];
s q0[0];
z q0[0];
s q0[0];
y q0[1];
s q0[1];
cx q0[0],q0[1];
h q0[0];
s q0[0];
sdg q0[1];
h q0[1];
cx q0[0],q0[1];
z q0[0];
h q0[0];
s q0[0];
h q0[1];
s q0[1];
x q0[1];
sdg q0[1];
h q0[1];
cx q0[0],q0[1];
h q0[0];
s q0[0];
s q0[0];
h q0[1];
s q0[1];
z q0[1];
sdg q0[1];
h q0[1];
cx q0[0],q0[1];
x q0[0];
s q0[0];
h q0[0];
s q0[0];
sdg q0[1];
h q0[1];
h q0[1];
cx q0[0],q0[1];
sdg q0[0];
h q0[0];
sdg q0[1];
h q0[1];
cx q0[0],q0[1];
h q0[0];
s q0[0];
z q0[0];
s q0[0];
x q0[1];
cx q0[0],q0[1];
h q0[0];
s q0[0];
x q0[0];
h q0[0];
h q0[1];
s q0[1];
z q0[1];
cx q0[0],q0[1];
sdg q0[0];
h q0[0];
h q0[1];
s q0[1];
cx q0[0],q0[1];
y q0[0];
h q0[0];
s q0[0];
h q0[1];
s q0[1];
cx q0[0],q0[1];
x q0[0];
h q0[0];
x q0[0];
s q0[0];
x q0[1];
s q0[1];
h q0[1];
s q0[1];
z q0[1];
cx q0[0],q0[1];
h q0[0];
s q0[0];
h q0[1];
s q0[1];
cx q0[0],q0[1];
s q0[0];
s q0[1];
cx q0[0],q0[1];
sdg q0[0];
h q0[0];
h q0[1];
s q0[1];
cx q0[0],q0[1];
h q0[0];
s q0[0];
x q0[0];
h q0[0];
s q0[0];
z q0[1];
s q0[1];
cx q0[0],q0[1];
sdg q0[0];
h q0[0];
h q0[1];
s q0[1];
cx q0[0],q0[1];
y q0[0];
sdg q0[0];
h q0[0];
y q0[1];
cx q0[0],q0[1];
sdg q0[0];
h q0[0];
sdg q0[1];
h q0[1];
cx q0[0],q0[1];
z q0[0];
sdg q0[0];
h q0[0];
h q0[1];
s q0[1];
z q0[1];
s q0[1];
h q0[1];
s q0[1];
cx q0[0],q0[1];
sdg q0[0];
h q0[0];
y q0[0];
sdg q0[0];
h q0[0];
z q0[1];
s q0[1];
h q0[1];
s q0[1];
cx q0[0],q0[1];
h q0[0];
s q0[0];
sdg q0[1];
h q0[1];
cx q0[0],q0[1];
z q0[0];
sdg q0[0];
h q0[0];
h q0[1];
s q0[1];
x q0[1];
s q0[1];
cx q0[0],q0[1];
sdg q0[0];
h q0[0];
h q0[1];
s q0[1];
cx q0[0],q0[1];
h q0[0];
s q0[0];
s q0[0];
h q0[0];
s q0[0];
y q0[1];
h q0[1];
cx q0[0],q0[1];
sdg q0[0];
h q0[0];
sdg q0[1];
h q0[1];
cx q0[0],q0[1];
z q0[0];
h q0[1];
s q0[1];
z q0[1];
h q0[1];
s q0[1];
cx q0[0],q0[1];
sdg q0[0];
h q0[0];
h q0[1];
s q0[1];
cx q0[0],q0[1];
y q0[0];
h q0[0];
x q0[1];
sdg q0[1];
h q0[1];
cx q0[0],q0[1];
sdg q0[0];
h q0[0];
sdg q0[1];
h q0[1];
cx q0[0],q0[1];
x q0[0];
h q0[1];
s q0[1];
z q0[1];
s q0[1];
cx q0[0],q0[1];
h q0[0];
s q0[0];
sdg q0[1];
h q0[1];
cx q0[0],q0[1];
y q0[0];
s q0[0];
h q0[1];
s q0[1];
y q0[1];
cx q0[0],q0[1];
sdg q0[0];
h q0[0];
h q0[1];
s q0[1];
cx q0[0],q0[1];
z q0[0];
h q0[0];
z q0[1];
cx q0[0],q0[1];
sdg q0[0];
h q0[0];
h q0[1];
s q0[1];
cx q0[0],q0[1];
h q0[0];
s q0[0];
h q0[0];
s q0[0];
y q0[1];
s q0[1];
cx q0[0],q0[1];
x q0[0];
h q0[0];
h q0[1];
s q0[1];
z q0[1];
h q0[1];
s q0[1];
cx q0[0],q0[1];
sdg q0[0];
h q0[0];
sdg q0[1];
h q0[1];
cx q0[0],q0[1];
h q0[0];
s q0[0];
z q0[0];
sdg q0[0];
h q0[0];
y q0[1];
cx q0[0],q0[1];
sdg q0[0];
h q0[0];
h q0[1];
s q0[1];
cx q0[0],q0[1];
z q0[0];
s q0[0];
z q0[1];
sdg q0[1];
h q0[1];
cx q0[0],q0[1];
h q0[0];
s q0[0];
z q0[0];
s q0[0];
x q0[1];
h q0[1];
s q0[1];
cx q0[0],q0[1];
sdg q0[0];
h q0[0];
h q0[1];
s q0[1];
cx q0[0],q0[1];
sdg q0[0];
h q0[0];
x q0[1];
h q0[1];
cx q0[0],q0[1];
sdg q0[0];
h q0[0];
h q0[1];
s q0[1];
cx q0[0],q0[1];
h q0[0];
s q0[0];
y q0[0];
sdg q0[0];
h q0[0];
x q0[1];
s q0[1];
cx q0[0],q0[1];
sdg q0[0];
h q0[0];
y q0[0];
h q0[0];
h q0[1];
s q0[1];
x q0[1];
sdg q0[1];
h q0[1];
cx q0[0],q0[1];
sdg q0[0];
h q0[0];
sdg q0[1];
h q0[1];
cx q0[0],q0[1];
h q0[0];
s q0[0];
x q0[0];
s q0[0];
h q0[0];
s q0[0];
z q0[1];
sdg q0[1];
h q0[1];
cx q0[0],q0[1];
sdg q0[0];
h q0[0];
sdg q0[1];
h q0[1];
cx q0[0],q0[1];
h q0[0];
s q0[0];
z q0[0];
s q0[0];
sdg q0[0];
h q0[0];
s q0[1];
h q0[1];
s q0[1];
s q0[1];
cx q0[0],q0[1];
sdg q0[0];
h q0[0];
x q0[1];
h q0[1];
cx q0[0],q0[1];
sdg q0[0];
h q0[0];
sdg q0[1];
h q0[1];
cx q0[0],q0[1];
h q0[0];
z q0[1];
s q0[1];
cx q0[0],q0[1];
x q0[0];
h q0[0];
z q0[0];
h q0[0];
y q0[1];
sdg q0[1];
h q0[1];
s q0[1];
h q0[1];
s q0[1];
cx q0[0],q0[1];
sdg q0[0];
h q0[0];
sdg q0[1];
h q0[1];
cx q0[0],q0[1];
h q0[1];
s q0[1];
x q0[1];
cx q0[0],q0[1];
sdg q0[0];
h q0[0];
sdg q0[1];
h q0[1];
cx q0[0],q0[1];
h q0[0];
s q0[0];
x q0[0];
cx q0[0],q0[1];
h q0[0];
s q0[0];
x q0[0];
h q0[0];
s q0[0];
z q0[1];
s q0[1];
h q0[1];
s q0[1];
cx q0[0],q0[1];
y q0[0];
h q0[0];
s q0[0];
sdg q0[1];
h q0[1];
z q0[1];
cx q0[0],q0[1];
sdg q0[0];
h q0[0];
sdg q0[1];
h q0[1];
cx q0[0],q0[1];
h q0[0];
s q0[0];
h q0[1];
s q0[1];
cx q0[0],q0[1];
x q0[0];
sdg q0[0];
h q0[0];
y q0[1];
sdg q0[1];
h q0[1];
cx q0[0],q0[1];
h q0[0];
s q0[0];
h q0[1];
s q0[1];
cx q0[0],q0[1];
y q0[0];
z q0[1];
h q0[1];
cx q0[0],q0[1];
sdg q0[0];
h q0[0];
x q0[0];
sdg q0[0];
h q0[0];
sdg q0[1];
h q0[1];
z q0[1];
s q0[1];
h q0[1];
s q0[1];
cx q0[0],q0[1];
sdg q0[0];
h q0[0];
sdg q0[1];
h q0[1];
cx q0[0],q0[1];
h q0[0];
s q0[0];
x q0[0];
s q0[0];
h q0[0];
s q0[0];
y q0[1];
h q0[1];
cx q0[0],q0[1];
h q0[0];
s q0[0];
sdg q0[1];
h q0[1];
cx q0[0],q0[1];
z q0[0];
s q0[0];
h q0[0];
s q0[0];
h q0[1];
s q0[1];
z q0[1];
h q0[1];
cx q0[0],q0[1];
h q0[0];
s q0[0];
x q0[0];
h q0[0];
s q0[0];
h q0[0];
h q0[1];
s q0[1];
y q0[1];
s q0[1];
z q0[1];
s q0[1];
h q0[1];
s q0[1];
cx q0[0],q0[1];
sdg q0[0];
h q0[0];
sdg q0[1];
h q0[1];
cx q0[0],q0[1];
x q0[0];
sdg q0[0];
h q0[0];
h q0[1];
s q0[1];
z q0[1];
cx q0[0],q0[1];
sdg q0[0];
h q0[0];
sdg q0[1];
h q0[1];
cx q0[0],q0[1];
h q0[0];
s q0[0];
x q0[0];
h q0[0];
h q0[1];
s q0[1];
x q0[1];
h q0[1];
cx q0[0],q0[1];
sdg q0[0];
h q0[0];
sdg q0[1];
h q0[1];
cx q0[0],q0[1];
h q0[0];
s q0[0];
x q0[0];
s q0[0];
h q0[0];
s q0[0];
h q0[1];
s q0[1];
y q0[1];
h q0[1];
s q0[1];
cx q0[0],q0[1];
h q0[0];
s q0[0];
sdg q0[1];
h q0[1];
cx q0[0],q0[1];
x q0[0];
sdg q0[0];
h q0[0];
y q0[1];
s q0[1];
h q0[1];
s q0[1];
cx q0[0],q0[1];
sdg q0[0];
h q0[0];
sdg q0[1];
h q0[1];
cx q0[0],q0[1];
h q0[0];
s q0[0];
h q0[1];
s q0[1];
cx q0[0],q0[1];
y q0[0];
h q0[0];
s q0[0];
x q0[1];
s q0[1];
cx q0[0],q0[1];
sdg q0[0];
h q0[0];
sdg q0[1];
h q0[1];
cx q0[0],q0[1];
z q0[0];
s q0[0];
h q0[0];
s q0[0];
h q0[1];
s q0[1];
x q0[1];
h q0[1];
s q0[1];
cx q0[0],q0[1];
sdg q0[0];
h q0[0];
h q0[1];
s q0[1];
cx q0[0],q0[1];
s q0[0];
y q0[1];
h q0[1];
s q0[1];
cx q0[0],q0[1];
h q0[0];
s q0[0];
x q0[0];
s q0[0];
s q0[1];
h q0[1];
s q0[1];
cx q0[0],q0[1];
sdg q0[0];
h q0[0];
h q0[1];
s q0[1];
cx q0[0],q0[1];
h q0[0];
s q0[0];
h q0[0];
s q0[0];
z q0[1];
sdg q0[1];
h q0[1];
cx q0[0],q0[1];
z q0[0];
h q0[0];
s q0[0];
h q0[1];
s q0[1];
z q0[1];
h q0[1];
cx q0[0],q0[1];
sdg q0[0];
h q0[0];
