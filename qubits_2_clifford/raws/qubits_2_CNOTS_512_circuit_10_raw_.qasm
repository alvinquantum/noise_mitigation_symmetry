OPENQASM 2.0;
include "qelib1.inc";
qreg q9[2];
y q9[0];
s q9[0];
z q9[1];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
z q9[0];
s q9[0];
s q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
y q9[0];
s q9[0];
h q9[0];
s q9[0];
h q9[1];
s q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
y q9[0];
s q9[0];
h q9[1];
s q9[1];
z q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
x q9[0];
s q9[0];
h q9[1];
s q9[1];
y q9[1];
s q9[1];
h q9[1];
s q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
y q9[0];
sdg q9[1];
h q9[1];
z q9[1];
s q9[1];
h q9[1];
s q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
y q9[0];
h q9[0];
s q9[0];
h q9[1];
s q9[1];
z q9[1];
h q9[1];
s q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
z q9[0];
sdg q9[0];
h q9[0];
h q9[1];
s q9[1];
z q9[1];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
h q9[1];
s q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
y q9[1];
h q9[1];
cx q9[0],q9[1];
y q9[0];
h q9[0];
s q9[0];
z q9[0];
s q9[0];
h q9[0];
s q9[0];
h q9[1];
s q9[1];
z q9[1];
z q9[1];
h q9[1];
s q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
z q9[1];
s q9[1];
cx q9[0],q9[1];
x q9[0];
s q9[0];
h q9[0];
s q9[0];
h q9[1];
s q9[1];
y q9[1];
s q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
h q9[1];
s q9[1];
y q9[1];
cx q9[0],q9[1];
y q9[0];
s q9[0];
h q9[0];
s q9[0];
z q9[1];
h q9[1];
s q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
x q9[0];
sdg q9[0];
h q9[0];
x q9[1];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
h q9[0];
s q9[0];
h q9[1];
s q9[1];
y q9[1];
h q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
y q9[0];
h q9[0];
x q9[1];
h q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
sdg q9[0];
h q9[0];
x q9[1];
h q9[1];
s q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
x q9[0];
s q9[0];
x q9[1];
h q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
z q9[0];
h q9[0];
h q9[1];
s q9[1];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
h q9[1];
s q9[1];
cx q9[0],q9[1];
x q9[0];
z q9[0];
s q9[0];
h q9[0];
s q9[0];
z q9[1];
y q9[1];
h q9[1];
cx q9[0],q9[1];
y q9[0];
h q9[0];
s q9[0];
sdg q9[1];
h q9[1];
h q9[1];
s q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
z q9[0];
h q9[0];
s q9[0];
h q9[1];
s q9[1];
x q9[1];
s q9[1];
h q9[1];
s q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
y q9[0];
s q9[0];
h q9[0];
s q9[0];
z q9[1];
h q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
z q9[0];
sdg q9[1];
h q9[1];
z q9[1];
s q9[1];
h q9[1];
s q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
h q9[1];
s q9[1];
cx q9[0],q9[1];
z q9[0];
sdg q9[0];
h q9[0];
y q9[1];
s q9[1];
h q9[1];
s q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
y q9[0];
s q9[0];
y q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
h q9[1];
s q9[1];
cx q9[0],q9[1];
y q9[0];
s q9[0];
h q9[0];
s q9[0];
z q9[0];
z q9[1];
s q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
z q9[0];
h q9[0];
s q9[0];
s q9[1];
h q9[1];
s q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
z q9[0];
sdg q9[0];
h q9[0];
h q9[1];
s q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
x q9[0];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
y q9[1];
h q9[1];
s q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
z q9[0];
s q9[0];
h q9[0];
s q9[0];
y q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
sdg q9[0];
h q9[0];
sdg q9[1];
h q9[1];
y q9[1];
h q9[1];
s q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
x q9[0];
h q9[0];
s q9[0];
z q9[1];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
x q9[0];
h q9[0];
sdg q9[1];
h q9[1];
y q9[1];
h q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
x q9[0];
h q9[0];
h q9[1];
s q9[1];
x q9[1];
h q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
h q9[1];
s q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
y q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
sdg q9[0];
h q9[0];
h q9[1];
s q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
z q9[0];
h q9[0];
sdg q9[1];
h q9[1];
x q9[1];
s q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
x q9[0];
h q9[0];
s q9[0];
x q9[1];
h q9[1];
s q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
z q9[0];
sdg q9[0];
h q9[0];
y q9[1];
s q9[1];
h q9[1];
s q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
x q9[0];
h q9[0];
s q9[0];
h q9[1];
s q9[1];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
s q9[0];
h q9[0];
s q9[0];
sdg q9[1];
h q9[1];
z q9[1];
h q9[1];
s q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
z q9[0];
h q9[0];
h q9[1];
s q9[1];
z q9[1];
s q9[1];
h q9[1];
s q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
h q9[1];
s q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
z q9[0];
h q9[0];
s q9[0];
x q9[1];
h q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
z q9[0];
s q9[0];
h q9[1];
s q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
h q9[1];
s q9[1];
cx q9[0],q9[1];
x q9[0];
s q9[0];
h q9[0];
s q9[0];
z q9[1];
h q9[1];
s q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
h q9[1];
s q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
z q9[0];
z q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
h q9[1];
s q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
x q9[0];
h q9[0];
s q9[0];
x q9[0];
sdg q9[0];
h q9[0];
y q9[1];
sdg q9[1];
h q9[1];
y q9[1];
s q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
y q9[0];
sdg q9[0];
h q9[0];
sdg q9[1];
h q9[1];
z q9[1];
h q9[1];
s q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
h q9[1];
s q9[1];
cx q9[0],q9[1];
s q9[0];
y q9[1];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
y q9[0];
s q9[0];
h q9[1];
s q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
h q9[1];
s q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
h q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
x q9[0];
h q9[0];
s q9[0];
sdg q9[1];
h q9[1];
h q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
x q9[0];
h q9[0];
h q9[1];
s q9[1];
z q9[1];
h q9[1];
s q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
s q9[0];
h q9[0];
s q9[0];
y q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
s q9[0];
h q9[0];
s q9[0];
h q9[1];
s q9[1];
z q9[1];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
x q9[0];
h q9[0];
s q9[0];
h q9[1];
s q9[1];
s q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
x q9[0];
h q9[1];
s q9[1];
z q9[1];
s q9[1];
h q9[1];
s q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
z q9[0];
h q9[0];
s q9[0];
h q9[1];
s q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
x q9[0];
sdg q9[0];
h q9[0];
h q9[1];
s q9[1];
z q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
y q9[0];
h q9[0];
s q9[0];
h q9[0];
s q9[0];
h q9[1];
s q9[1];
y q9[1];
h q9[1];
y q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
y q9[0];
s q9[0];
h q9[0];
s q9[0];
h q9[1];
s q9[1];
y q9[1];
h q9[1];
s q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
x q9[0];
h q9[1];
s q9[1];
y q9[1];
s q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
h q9[1];
s q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
x q9[0];
s q9[0];
y q9[1];
s q9[1];
h q9[1];
s q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
x q9[0];
h q9[0];
h q9[1];
s q9[1];
y q9[1];
h q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
z q9[0];
s q9[0];
h q9[0];
s q9[0];
h q9[1];
s q9[1];
x q9[1];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
y q9[0];
s q9[0];
h q9[0];
s q9[0];
x q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
h q9[1];
s q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
z q9[0];
x q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
x q9[0];
sdg q9[0];
h q9[0];
h q9[1];
s q9[1];
x q9[1];
h q9[1];
s q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
y q9[0];
s q9[0];
sdg q9[1];
h q9[1];
x q9[1];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
z q9[0];
h q9[1];
s q9[1];
z q9[1];
s q9[1];
h q9[1];
s q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
s q9[0];
h q9[0];
s q9[0];
h q9[1];
s q9[1];
h q9[1];
s q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
sdg q9[0];
h q9[0];
x q9[1];
s q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
sdg q9[1];
h q9[1];
x q9[1];
cx q9[0],q9[1];
s q9[0];
h q9[0];
s q9[0];
h q9[1];
s q9[1];
z q9[1];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
y q9[0];
s q9[0];
sdg q9[1];
h q9[1];
h q9[1];
s q9[1];
cx q9[0],q9[1];
s q9[0];
sdg q9[1];
h q9[1];
s q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
h q9[1];
s q9[1];
cx q9[0],q9[1];
x q9[0];
sdg q9[0];
h q9[0];
z q9[1];
h q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
y q9[0];
s q9[0];
h q9[0];
s q9[0];
y q9[1];
s q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
y q9[0];
h q9[0];
s q9[0];
sdg q9[1];
h q9[1];
x q9[1];
s q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
y q9[0];
h q9[0];
z q9[0];
h q9[1];
s q9[1];
s q9[1];
h q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
x q9[0];
s q9[0];
z q9[1];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
sdg q9[0];
h q9[0];
sdg q9[1];
h q9[1];
s q9[1];
h q9[1];
s q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
h q9[1];
s q9[1];
cx q9[0],q9[1];
x q9[0];
h q9[0];
s q9[0];
y q9[1];
h q9[1];
s q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
z q9[0];
x q9[1];
h q9[1];
s q9[1];
cx q9[0],q9[1];
x q9[0];
h q9[0];
x q9[1];
s q9[1];
h q9[1];
s q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
h q9[1];
s q9[1];
cx q9[0],q9[1];
x q9[0];
s q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
x q9[0];
s q9[0];
h q9[0];
s q9[0];
z q9[0];
h q9[0];
s q9[0];
sdg q9[1];
h q9[1];
z q9[1];
sdg q9[1];
h q9[1];
x q9[1];
h q9[1];
s q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
z q9[0];
x q9[1];
s q9[1];
h q9[1];
s q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
h q9[0];
h q9[1];
s q9[1];
x q9[1];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
h q9[1];
s q9[1];
cx q9[0],q9[1];
z q9[0];
h q9[0];
y q9[1];
s q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
y q9[0];
s q9[0];
h q9[0];
s q9[0];
sdg q9[1];
h q9[1];
z q9[1];
s q9[1];
h q9[1];
s q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
h q9[1];
s q9[1];
cx q9[0],q9[1];
z q9[0];
sdg q9[0];
h q9[0];
z q9[1];
s q9[1];
h q9[1];
s q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
h q9[1];
s q9[1];
cx q9[0],q9[1];
z q9[0];
s q9[0];
h q9[0];
s q9[0];
z q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
h q9[0];
s q9[0];
sdg q9[1];
h q9[1];
y q9[1];
h q9[1];
s q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
h q9[1];
s q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
sdg q9[0];
h q9[0];
y q9[1];
s q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
h q9[1];
s q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
z q9[0];
h q9[0];
x q9[1];
h q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
z q9[1];
h q9[1];
s q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
z q9[0];
y q9[1];
s q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
h q9[1];
s q9[1];
cx q9[0],q9[1];
z q9[0];
z q9[1];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
sdg q9[0];
h q9[0];
h q9[1];
s q9[1];
z q9[1];
h q9[1];
s q9[1];
cx q9[0],q9[1];
x q9[0];
s q9[0];
x q9[1];
s q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
h q9[1];
s q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
sdg q9[0];
h q9[0];
x q9[1];
h q9[1];
s q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
x q9[0];
sdg q9[0];
h q9[0];
h q9[1];
s q9[1];
h q9[1];
s q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
z q9[0];
h q9[1];
s q9[1];
s q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
s q9[0];
s q9[1];
h q9[1];
s q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
h q9[1];
s q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
z q9[0];
sdg q9[0];
h q9[0];
y q9[1];
h q9[1];
s q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
y q9[0];
h q9[0];
s q9[0];
y q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
x q9[0];
s q9[0];
sdg q9[1];
h q9[1];
s q9[1];
h q9[1];
s q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
x q9[0];
h q9[0];
s q9[0];
h q9[1];
s q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
x q9[0];
s q9[0];
h q9[0];
s q9[0];
sdg q9[1];
h q9[1];
h q9[1];
s q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
h q9[1];
s q9[1];
cx q9[0],q9[1];
z q9[0];
sdg q9[0];
h q9[0];
s q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
z q9[0];
sdg q9[0];
h q9[0];
x q9[1];
h q9[1];
s q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
s q9[0];
z q9[1];
s q9[1];
h q9[1];
s q9[1];
cx q9[0],q9[1];
y q9[0];
h q9[1];
s q9[1];
z q9[1];
s q9[1];
h q9[1];
s q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
y q9[0];
x q9[1];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
h q9[1];
s q9[1];
cx q9[0],q9[1];
x q9[0];
h q9[1];
s q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
h q9[1];
s q9[1];
cx q9[0],q9[1];
x q9[1];
s q9[1];
h q9[1];
s q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
s q9[1];
cx q9[0],q9[1];
x q9[0];
h q9[0];
z q9[1];
s q9[1];
h q9[1];
s q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
h q9[1];
s q9[1];
cx q9[0],q9[1];
y q9[0];
y q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
h q9[1];
s q9[1];
cx q9[0],q9[1];
x q9[0];
h q9[0];
s q9[0];
h q9[1];
s q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
y q9[0];
h q9[0];
s q9[0];
s q9[0];
h q9[0];
s q9[0];
z q9[1];
s q9[1];
x q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
sdg q9[0];
h q9[0];
h q9[1];
s q9[1];
cx q9[0],q9[1];
z q9[0];
s q9[0];
h q9[0];
s q9[0];
y q9[1];
h q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
z q9[0];
s q9[0];
h q9[0];
s q9[0];
h q9[1];
s q9[1];
y q9[1];
h q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
h q9[1];
s q9[1];
cx q9[0],q9[1];
y q9[0];
sdg q9[0];
h q9[0];
z q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
z q9[0];
s q9[0];
h q9[0];
s q9[0];
h q9[1];
s q9[1];
y q9[1];
h q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
h q9[1];
s q9[1];
cx q9[0],q9[1];
z q9[1];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
x q9[0];
sdg q9[0];
h q9[0];
y q9[0];
sdg q9[0];
h q9[0];
h q9[1];
s q9[1];
x q9[1];
h q9[1];
s q9[1];
x q9[1];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
y q9[0];
h q9[0];
s q9[0];
y q9[1];
s q9[1];
h q9[1];
s q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
z q9[0];
h q9[0];
y q9[1];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
h q9[0];
x q9[1];
s q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
x q9[0];
h q9[0];
s q9[0];
h q9[1];
s q9[1];
y q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
s q9[0];
h q9[1];
s q9[1];
x q9[1];
h q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
h q9[1];
s q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
z q9[1];
h q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
h q9[1];
s q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
y q9[0];
sdg q9[0];
h q9[0];
y q9[1];
h q9[1];
s q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
y q9[0];
sdg q9[0];
h q9[0];
z q9[1];
s q9[1];
h q9[1];
s q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
z q9[1];
s q9[1];
h q9[1];
s q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
x q9[0];
sdg q9[0];
h q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
sdg q9[0];
h q9[0];
h q9[1];
s q9[1];
x q9[1];
s q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
x q9[0];
y q9[1];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
h q9[1];
s q9[1];
cx q9[0],q9[1];
z q9[0];
s q9[0];
h q9[0];
s q9[0];
x q9[1];
s q9[1];
h q9[1];
s q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
s q9[0];
h q9[0];
s q9[0];
h q9[1];
s q9[1];
y q9[1];
s q9[1];
h q9[1];
s q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
y q9[0];
sdg q9[0];
h q9[0];
y q9[1];
h q9[1];
s q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
z q9[0];
h q9[0];
s q9[0];
h q9[1];
s q9[1];
y q9[1];
h q9[1];
cx q9[0],q9[1];
x q9[0];
s q9[0];
h q9[0];
s q9[0];
h q9[1];
s q9[1];
x q9[1];
h q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
z q9[0];
h q9[1];
s q9[1];
z q9[1];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
y q9[0];
sdg q9[0];
h q9[0];
h q9[1];
s q9[1];
x q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
y q9[0];
h q9[0];
s q9[0];
sdg q9[1];
h q9[1];
s q9[1];
h q9[1];
s q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
z q9[0];
h q9[1];
s q9[1];
z q9[1];
h q9[1];
s q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
z q9[0];
h q9[0];
s q9[0];
y q9[1];
s q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
y q9[0];
h q9[0];
s q9[0];
sdg q9[1];
h q9[1];
s q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
z q9[0];
h q9[0];
s q9[0];
h q9[1];
s q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
z q9[0];
s q9[0];
h q9[0];
s q9[0];
h q9[1];
s q9[1];
y q9[1];
s q9[1];
h q9[1];
s q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
h q9[1];
s q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
z q9[0];
sdg q9[0];
h q9[0];
cx q9[0],q9[1];
h q9[0];
s q9[0];
z q9[0];
h q9[0];
h q9[1];
s q9[1];
z q9[1];
s q9[1];
h q9[1];
s q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
s q9[0];
h q9[1];
s q9[1];
s q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
z q9[0];
s q9[0];
z q9[1];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
h q9[1];
s q9[1];
y q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
h q9[1];
s q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
z q9[1];
h q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
s q9[0];
sdg q9[1];
h q9[1];
x q9[1];
cx q9[0],q9[1];
y q9[0];
sdg q9[1];
h q9[1];
y q9[1];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
y q9[0];
s q9[0];
sdg q9[1];
h q9[1];
s q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
y q9[0];
h q9[0];
s q9[0];
h q9[1];
s q9[1];
y q9[1];
h q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
h q9[1];
s q9[1];
cx q9[0],q9[1];
x q9[0];
s q9[0];
z q9[1];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
z q9[0];
h q9[1];
s q9[1];
y q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
h q9[1];
s q9[1];
cx q9[0],q9[1];
x q9[0];
s q9[0];
x q9[1];
s q9[1];
h q9[1];
s q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
sdg q9[1];
h q9[1];
y q9[1];
h q9[1];
s q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
h q9[1];
s q9[1];
cx q9[0],q9[1];
z q9[0];
h q9[0];
x q9[0];
z q9[1];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
s q9[0];
z q9[1];
s q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
z q9[0];
h q9[0];
s q9[0];
sdg q9[1];
h q9[1];
z q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
h q9[1];
s q9[1];
cx q9[0],q9[1];
x q9[0];
h q9[0];
y q9[1];
s q9[1];
h q9[1];
s q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
h q9[1];
s q9[1];
cx q9[0],q9[1];
z q9[0];
s q9[0];
x q9[0];
h q9[0];
s q9[0];
x q9[1];
s q9[1];
x q9[1];
h q9[1];
cx q9[0],q9[1];
x q9[0];
sdg q9[1];
h q9[1];
y q9[1];
h q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
z q9[0];
h q9[0];
y q9[0];
s q9[0];
h q9[0];
s q9[0];
x q9[1];
z q9[1];
s q9[1];
h q9[1];
s q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
x q9[0];
sdg q9[0];
h q9[0];
h q9[1];
s q9[1];
x q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
h q9[1];
s q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
x q9[0];
s q9[0];
h q9[0];
s q9[0];
x q9[1];
s q9[1];
h q9[1];
s q9[1];
cx q9[0],q9[1];
z q9[0];
sdg q9[0];
h q9[0];
h q9[1];
s q9[1];
x q9[1];
s q9[1];
h q9[1];
s q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
y q9[0];
h q9[0];
s q9[0];
x q9[1];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
z q9[0];
s q9[0];
h q9[0];
s q9[0];
sdg q9[0];
h q9[0];
sdg q9[1];
h q9[1];
x q9[1];
s q9[1];
x q9[1];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
y q9[0];
h q9[0];
s q9[0];
h q9[1];
s q9[1];
s q9[1];
h q9[1];
s q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
sdg q9[0];
h q9[0];
h q9[1];
s q9[1];
y q9[1];
h q9[1];
s q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
y q9[0];
h q9[0];
s q9[0];
h q9[1];
s q9[1];
x q9[1];
h q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
h q9[1];
s q9[1];
cx q9[0],q9[1];
x q9[0];
h q9[0];
s q9[0];
y q9[1];
s q9[1];
h q9[1];
s q9[1];
cx q9[0],q9[1];
z q9[0];
h q9[0];
s q9[0];
h q9[1];
s q9[1];
x q9[1];
h q9[1];
s q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
h q9[0];
h q9[1];
s q9[1];
h q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
h q9[1];
s q9[1];
x q9[1];
s q9[1];
h q9[1];
s q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
h q9[1];
s q9[1];
cx q9[0],q9[1];
z q9[0];
h q9[1];
s q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
s q9[0];
h q9[1];
s q9[1];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
y q9[0];
h q9[0];
s q9[0];
h q9[1];
s q9[1];
z q9[1];
h q9[1];
s q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
x q9[0];
s q9[0];
h q9[0];
s q9[0];
y q9[1];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
z q9[0];
z q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
sdg q9[0];
h q9[0];
x q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
x q9[0];
h q9[0];
s q9[0];
h q9[1];
s q9[1];
h q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
x q9[0];
h q9[0];
sdg q9[1];
h q9[1];
h q9[1];
s q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
y q9[0];
s q9[0];
h q9[0];
s q9[0];
h q9[1];
s q9[1];
z q9[1];
s q9[1];
h q9[1];
s q9[1];
cx q9[0],q9[1];
x q9[0];
h q9[0];
z q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
h q9[1];
s q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
y q9[1];
h q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
h q9[0];
s q9[0];
y q9[1];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
y q9[0];
s q9[0];
sdg q9[1];
h q9[1];
x q9[1];
s q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
y q9[0];
h q9[1];
s q9[1];
z q9[1];
s q9[1];
h q9[1];
s q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
y q9[0];
s q9[0];
h q9[1];
s q9[1];
h q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
y q9[0];
s q9[0];
h q9[0];
s q9[0];
h q9[1];
s q9[1];
x q9[1];
s q9[1];
h q9[1];
s q9[1];
cx q9[0],q9[1];
y q9[0];
sdg q9[0];
h q9[0];
z q9[0];
h q9[1];
s q9[1];
z q9[1];
y q9[1];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
y q9[0];
s q9[0];
h q9[0];
s q9[0];
h q9[1];
s q9[1];
h q9[1];
s q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
x q9[0];
sdg q9[0];
h q9[0];
x q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
h q9[1];
s q9[1];
cx q9[0],q9[1];
y q9[0];
s q9[0];
h q9[0];
s q9[0];
s q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
h q9[1];
s q9[1];
cx q9[0],q9[1];
z q9[0];
h q9[0];
s q9[0];
z q9[1];
h q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
sdg q9[0];
h q9[0];
z q9[1];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
h q9[1];
s q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
z q9[0];
s q9[0];
y q9[1];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
h q9[1];
s q9[1];
cx q9[0],q9[1];
s q9[1];
h q9[1];
s q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
x q9[0];
h q9[0];
y q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
h q9[1];
s q9[1];
cx q9[0],q9[1];
y q9[0];
s q9[0];
h q9[0];
s q9[0];
z q9[1];
s q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
s q9[0];
h q9[0];
s q9[0];
h q9[1];
cx q9[0],q9[1];
x q9[0];
sdg q9[0];
h q9[0];
y q9[0];
s q9[0];
h q9[0];
s q9[0];
x q9[1];
s q9[1];
h q9[1];
s q9[1];
x q9[1];
h q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
x q9[0];
sdg q9[0];
h q9[0];
y q9[1];
h q9[1];
s q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
h q9[1];
s q9[1];
cx q9[0],q9[1];
y q9[0];
h q9[0];
s q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
x q9[0];
s q9[0];
h q9[0];
s q9[0];
sdg q9[1];
h q9[1];
x q9[1];
s q9[1];
h q9[1];
s q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
y q9[0];
h q9[0];
h q9[1];
s q9[1];
x q9[1];
s q9[1];
cx q9[0],q9[1];
z q9[0];
h q9[1];
s q9[1];
z q9[1];
h q9[1];
s q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
h q9[1];
s q9[1];
cx q9[0],q9[1];
y q9[0];
h q9[1];
s q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
z q9[0];
h q9[0];
s q9[0];
h q9[1];
s q9[1];
h q9[1];
s q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
z q9[0];
sdg q9[0];
h q9[0];
y q9[1];
h q9[1];
s q9[1];
cx q9[0],q9[1];
y q9[0];
h q9[0];
h q9[1];
s q9[1];
z q9[1];
h q9[1];
s q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
x q9[0];
s q9[0];
sdg q9[1];
h q9[1];
x q9[1];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
z q9[0];
z q9[1];
s q9[1];
h q9[1];
s q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
y q9[0];
s q9[0];
h q9[1];
s q9[1];
x q9[1];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
y q9[0];
h q9[0];
h q9[1];
s q9[1];
z q9[1];
h q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
z q9[0];
s q9[0];
h q9[0];
s q9[0];
x q9[1];
s q9[1];
h q9[1];
s q9[1];
cx q9[0],q9[1];
y q9[1];
s q9[1];
cx q9[0],q9[1];
y q9[0];
sdg q9[0];
h q9[0];
x q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
z q9[0];
sdg q9[0];
h q9[0];
h q9[1];
s q9[1];
x q9[1];
cx q9[0],q9[1];
s q9[0];
h q9[0];
s q9[0];
sdg q9[1];
h q9[1];
x q9[1];
s q9[1];
h q9[1];
s q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
y q9[0];
h q9[0];
h q9[1];
s q9[1];
x q9[1];
h q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
y q9[0];
s q9[0];
h q9[0];
s q9[0];
h q9[1];
s q9[1];
z q9[1];
h q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
y q9[0];
s q9[0];
h q9[1];
s q9[1];
z q9[1];
s q9[1];
h q9[1];
s q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
z q9[0];
s q9[0];
z q9[1];
h q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
x q9[0];
sdg q9[0];
h q9[0];
h q9[1];
s q9[1];
y q9[1];
s q9[1];
h q9[1];
s q9[1];
cx q9[0],q9[1];
z q9[0];
sdg q9[0];
h q9[0];
y q9[1];
h q9[1];
s q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
x q9[0];
s q9[0];
h q9[1];
s q9[1];
s q9[1];
h q9[1];
s q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
y q9[0];
sdg q9[0];
h q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
h q9[1];
s q9[1];
cx q9[0],q9[1];
s q9[0];
h q9[0];
s q9[0];
x q9[1];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
z q9[0];
h q9[1];
s q9[1];
x q9[1];
cx q9[0],q9[1];
y q9[0];
h q9[0];
sdg q9[1];
h q9[1];
h q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
h q9[1];
s q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
h q9[1];
s q9[1];
cx q9[0],q9[1];
x q9[0];
h q9[0];
h q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
z q9[0];
s q9[0];
x q9[1];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
sdg q9[0];
h q9[0];
h q9[1];
s q9[1];
y q9[1];
s q9[1];
h q9[1];
s q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
h q9[1];
s q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
y q9[0];
s q9[0];
h q9[0];
s q9[0];
x q9[1];
h q9[1];
s q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
x q9[0];
h q9[0];
s q9[0];
h q9[1];
s q9[1];
z q9[1];
s q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
x q9[0];
s q9[0];
h q9[0];
s q9[0];
y q9[1];
h q9[1];
z q9[1];
cx q9[0],q9[1];
x q9[0];
s q9[0];
h q9[0];
s q9[0];
cx q9[0],q9[1];
h q9[0];
s q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
x q9[0];
h q9[0];
s q9[0];
h q9[1];
s q9[1];
z q9[1];
h q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
y q9[0];
h q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
h q9[1];
s q9[1];
cx q9[0],q9[1];
z q9[0];
sdg q9[0];
h q9[0];
x q9[0];
s q9[0];
h q9[1];
s q9[1];
x q9[1];
h q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
x q9[0];
s q9[0];
z q9[1];
s q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
z q9[0];
sdg q9[0];
h q9[0];
sdg q9[1];
h q9[1];
z q9[1];
h q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
x q9[0];
s q9[0];
y q9[1];
s q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
h q9[1];
s q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
x q9[0];
sdg q9[0];
h q9[0];
z q9[1];
h q9[1];
s q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
x q9[0];
s q9[0];
h q9[0];
s q9[0];
h q9[1];
s q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
y q9[0];
h q9[0];
s q9[0];
h q9[1];
s q9[1];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
x q9[0];
h q9[1];
s q9[1];
y q9[1];
s q9[1];
h q9[1];
s q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
z q9[0];
s q9[0];
y q9[1];
h q9[1];
s q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
h q9[1];
s q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
x q9[0];
s q9[0];
h q9[0];
s q9[0];
h q9[1];
s q9[1];
y q9[1];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
z q9[0];
sdg q9[0];
h q9[0];
x q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
h q9[1];
s q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
h q9[1];
s q9[1];
cx q9[0],q9[1];
y q9[1];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
x q9[0];
s q9[0];
h q9[1];
s q9[1];
y q9[1];
h q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
h q9[1];
s q9[1];
cx q9[0],q9[1];
y q9[0];
h q9[0];
x q9[1];
h q9[1];
s q9[1];
cx q9[0],q9[1];
s q9[0];
h q9[0];
s q9[0];
h q9[1];
s q9[1];
z q9[1];
s q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
x q9[0];
h q9[1];
s q9[1];
x q9[1];
s q9[1];
cx q9[0],q9[1];
sdg q9[1];
h q9[1];
y q9[1];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
y q9[0];
h q9[0];
h q9[1];
s q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
y q9[0];
h q9[0];
s q9[0];
y q9[1];
h q9[1];
s q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
y q9[0];
sdg q9[0];
h q9[0];
sdg q9[1];
h q9[1];
h q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
h q9[1];
s q9[1];
cx q9[0],q9[1];
x q9[0];
z q9[1];
h q9[1];
s q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
h q9[0];
h q9[1];
s q9[1];
y q9[1];
h q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
y q9[0];
h q9[0];
z q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
y q9[0];
x q9[1];
h q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
sdg q9[0];
h q9[0];
cx q9[0],q9[1];
h q9[0];
s q9[0];
z q9[0];
s q9[0];
h q9[0];
s q9[0];
h q9[1];
s q9[1];
x q9[1];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
h q9[0];
h q9[1];
s q9[1];
x q9[1];
h q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
y q9[0];
h q9[0];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
h q9[1];
s q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
cx q9[0],q9[1];
h q9[0];
s q9[0];
y q9[0];
s q9[0];
h q9[0];
s q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
y q9[0];
s q9[0];
sdg q9[1];
h q9[1];
y q9[1];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
z q9[0];
h q9[0];
s q9[0];
h q9[1];
s q9[1];
z q9[1];
h q9[1];
s q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
sdg q9[0];
h q9[0];
h q9[1];
s q9[1];
z q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
h q9[0];
h q9[1];
s q9[1];
z q9[1];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
y q9[0];
s q9[0];
s q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
z q9[0];
sdg q9[0];
h q9[0];
h q9[1];
s q9[1];
y q9[1];
s q9[1];
h q9[1];
s q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
y q9[0];
s q9[0];
h q9[0];
s q9[0];
sdg q9[1];
h q9[1];
z q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
y q9[0];
sdg q9[0];
h q9[0];
x q9[1];
h q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
h q9[1];
s q9[1];
y q9[1];
s q9[1];
cx q9[0],q9[1];
x q9[0];
h q9[0];
sdg q9[1];
h q9[1];
y q9[1];
s q9[1];
h q9[1];
s q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
y q9[0];
h q9[1];
s q9[1];
x q9[1];
h q9[1];
s q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
h q9[1];
s q9[1];
cx q9[0],q9[1];
s q9[0];
h q9[0];
s q9[0];
cx q9[0],q9[1];
x q9[0];
sdg q9[0];
h q9[0];
h q9[1];
s q9[1];
s q9[1];
h q9[1];
s q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
x q9[0];
sdg q9[1];
h q9[1];
y q9[1];
h q9[1];
s q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
x q9[0];
z q9[1];
s q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
y q9[0];
h q9[0];
z q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
x q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
sdg q9[0];
h q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
h q9[1];
s q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
y q9[0];
cx q9[0],q9[1];
h q9[0];
s q9[0];
h q9[1];
s q9[1];
cx q9[0],q9[1];
x q9[0];
h q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
x q9[0];
s q9[0];
h q9[0];
s q9[0];
h q9[1];
s q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
h q9[0];
s q9[0];
s q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
h q9[1];
s q9[1];
cx q9[0],q9[1];
z q9[0];
sdg q9[0];
h q9[0];
x q9[1];
s q9[1];
h q9[1];
s q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
z q9[0];
h q9[0];
s q9[0];
z q9[1];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
z q9[0];
sdg q9[0];
h q9[0];
sdg q9[1];
h q9[1];
z q9[1];
h q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
h q9[1];
s q9[1];
x q9[1];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
h q9[0];
s q9[0];
z q9[1];
s q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
h q9[0];
h q9[1];
s q9[1];
z q9[1];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
h q9[1];
s q9[1];
x q9[1];
h q9[1];
s q9[1];
cx q9[0],q9[1];
z q9[0];
s q9[0];
h q9[0];
s q9[0];
sdg q9[1];
h q9[1];
x q9[1];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
z q9[0];
h q9[0];
h q9[1];
s q9[1];
z q9[1];
h q9[1];
s q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
h q9[0];
x q9[1];
cx q9[0],q9[1];
y q9[0];
sdg q9[0];
h q9[0];
h q9[1];
s q9[1];
z q9[1];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
y q9[0];
s q9[0];
h q9[0];
s q9[0];
x q9[1];
h q9[1];
s q9[1];
cx q9[0],q9[1];
y q9[0];
y q9[1];
h q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
x q9[0];
h q9[0];
x q9[1];
s q9[1];
h q9[1];
s q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
z q9[0];
h q9[0];
h q9[1];
s q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
y q9[0];
sdg q9[0];
h q9[0];
h q9[1];
s q9[1];
y q9[1];
s q9[1];
h q9[1];
s q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
y q9[0];
s q9[0];
h q9[0];
s q9[0];
sdg q9[1];
h q9[1];
s q9[1];
h q9[1];
s q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
sdg q9[1];
h q9[1];
cx q9[0],q9[1];
h q9[0];
s q9[0];
x q9[0];
s q9[0];
h q9[1];
s q9[1];
y q9[1];
s q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
x q9[0];
h q9[0];
sdg q9[1];
h q9[1];
s q9[1];
cx q9[0],q9[1];
sdg q9[0];
h q9[0];
sdg q9[1];
h q9[1];
