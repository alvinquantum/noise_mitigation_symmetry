OPENQASM 2.0;
include "qelib1.inc";
qreg q16[2];
y q16[0];
h q16[0];
h q16[1];
cx q16[0],q16[1];
y q16[0];
h q16[0];
y q16[0];
sdg q16[0];
h q16[0];
x q16[1];
s q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
h q16[0];
h q16[1];
s q16[1];
y q16[1];
s q16[1];
h q16[1];
s q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
y q16[0];
h q16[0];
h q16[1];
s q16[1];
z q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
h q16[1];
s q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
z q16[0];
sdg q16[0];
h q16[0];
z q16[1];
cx q16[0],q16[1];
z q16[0];
sdg q16[0];
h q16[0];
h q16[1];
s q16[1];
x q16[1];
h q16[1];
s q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
s q16[0];
x q16[1];
h q16[1];
s q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
y q16[0];
s q16[0];
h q16[0];
s q16[0];
x q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
h q16[1];
s q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
s q16[0];
x q16[1];
s q16[1];
h q16[1];
s q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
z q16[0];
h q16[1];
s q16[1];
y q16[1];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
z q16[0];
s q16[0];
h q16[1];
s q16[1];
x q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
h q16[0];
h q16[1];
s q16[1];
y q16[1];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
x q16[0];
sdg q16[0];
h q16[0];
h q16[1];
s q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
h q16[1];
s q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
x q16[0];
s q16[0];
h q16[0];
s q16[0];
y q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
z q16[0];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
x q16[0];
s q16[0];
h q16[0];
s q16[0];
y q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
h q16[1];
s q16[1];
cx q16[0],q16[1];
z q16[0];
h q16[0];
y q16[1];
h q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
h q16[1];
s q16[1];
cx q16[0],q16[1];
y q16[0];
s q16[0];
h q16[0];
s q16[0];
s q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
x q16[0];
h q16[0];
s q16[0];
x q16[0];
h q16[0];
z q16[0];
s q16[0];
h q16[0];
s q16[0];
h q16[1];
s q16[1];
h q16[1];
s q16[1];
z q16[1];
h q16[1];
s q16[1];
y q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
h q16[1];
s q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
x q16[1];
s q16[1];
h q16[1];
s q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
z q16[0];
s q16[0];
y q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
h q16[1];
s q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
x q16[0];
sdg q16[0];
h q16[0];
z q16[1];
s q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
z q16[0];
h q16[0];
z q16[1];
h q16[1];
cx q16[0],q16[1];
z q16[0];
h q16[0];
s q16[0];
y q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
h q16[1];
s q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
y q16[0];
x q16[1];
s q16[1];
h q16[1];
s q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
h q16[1];
s q16[1];
cx q16[0],q16[1];
z q16[0];
h q16[0];
y q16[1];
h q16[1];
s q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
x q16[0];
s q16[0];
h q16[0];
s q16[0];
z q16[1];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
z q16[0];
h q16[0];
s q16[0];
sdg q16[1];
h q16[1];
x q16[1];
s q16[1];
h q16[1];
s q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
y q16[0];
h q16[0];
s q16[0];
x q16[1];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
y q16[0];
h q16[0];
h q16[1];
s q16[1];
y q16[1];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
h q16[1];
s q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
x q16[0];
sdg q16[0];
h q16[0];
y q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
x q16[0];
s q16[0];
s q16[0];
h q16[0];
s q16[0];
h q16[1];
s q16[1];
y q16[1];
s q16[1];
h q16[1];
s q16[1];
x q16[1];
h q16[1];
s q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
h q16[1];
s q16[1];
cx q16[0],q16[1];
x q16[0];
s q16[0];
h q16[0];
s q16[0];
y q16[1];
h q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
x q16[0];
sdg q16[0];
h q16[0];
h q16[1];
s q16[1];
z q16[1];
h q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
s q16[0];
h q16[0];
s q16[0];
sdg q16[1];
h q16[1];
y q16[1];
h q16[1];
s q16[1];
cx q16[0],q16[1];
x q16[0];
z q16[1];
h q16[1];
cx q16[0],q16[1];
s q16[0];
z q16[1];
h q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
h q16[0];
s q16[0];
h q16[1];
s q16[1];
z q16[1];
h q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
z q16[0];
h q16[1];
s q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
x q16[0];
s q16[0];
h q16[0];
s q16[0];
h q16[1];
s q16[1];
s q16[1];
h q16[1];
s q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
h q16[1];
s q16[1];
cx q16[0],q16[1];
y q16[0];
s q16[0];
x q16[1];
s q16[1];
h q16[1];
s q16[1];
cx q16[0],q16[1];
z q16[0];
h q16[0];
s q16[0];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
y q16[0];
s q16[0];
h q16[0];
s q16[0];
h q16[1];
s q16[1];
y q16[1];
h q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
h q16[1];
s q16[1];
cx q16[0],q16[1];
x q16[0];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
s q16[0];
h q16[0];
s q16[0];
x q16[1];
h q16[1];
s q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
h q16[1];
s q16[1];
cx q16[0],q16[1];
z q16[0];
h q16[0];
s q16[0];
y q16[1];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
x q16[0];
s q16[0];
sdg q16[1];
h q16[1];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
h q16[1];
s q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
y q16[0];
sdg q16[0];
h q16[0];
y q16[1];
s q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
y q16[0];
s q16[0];
h q16[0];
s q16[0];
h q16[1];
s q16[1];
h q16[1];
s q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
z q16[0];
s q16[0];
h q16[1];
s q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
h q16[1];
s q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
y q16[0];
s q16[0];
z q16[1];
s q16[1];
h q16[1];
s q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
y q16[0];
h q16[0];
h q16[1];
s q16[1];
s q16[1];
h q16[1];
s q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
s q16[0];
h q16[0];
s q16[0];
sdg q16[1];
h q16[1];
s q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
h q16[1];
s q16[1];
cx q16[0],q16[1];
x q16[0];
h q16[0];
s q16[0];
x q16[1];
s q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
sdg q16[0];
h q16[0];
h q16[1];
s q16[1];
x q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
h q16[1];
s q16[1];
x q16[1];
h q16[1];
s q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
z q16[0];
s q16[0];
h q16[1];
s q16[1];
x q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
z q16[0];
s q16[0];
h q16[0];
s q16[0];
z q16[1];
s q16[1];
h q16[1];
s q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
h q16[1];
s q16[1];
cx q16[0],q16[1];
z q16[0];
s q16[0];
y q16[1];
s q16[1];
h q16[1];
s q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
y q16[0];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
z q16[0];
s q16[0];
h q16[1];
s q16[1];
x q16[1];
h q16[1];
s q16[1];
cx q16[0],q16[1];
x q16[0];
s q16[0];
h q16[0];
s q16[0];
y q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
y q16[0];
s q16[0];
h q16[0];
s q16[0];
h q16[1];
s q16[1];
y q16[1];
h q16[1];
x q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
h q16[1];
s q16[1];
cx q16[0],q16[1];
y q16[0];
sdg q16[0];
h q16[0];
z q16[1];
s q16[1];
h q16[1];
s q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
s q16[0];
s q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
x q16[0];
s q16[0];
h q16[1];
s q16[1];
h q16[1];
cx q16[0],q16[1];
x q16[0];
sdg q16[0];
h q16[0];
h q16[1];
s q16[1];
y q16[1];
s q16[1];
h q16[1];
s q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
y q16[0];
h q16[0];
s q16[0];
h q16[1];
s q16[1];
y q16[1];
s q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
x q16[0];
s q16[0];
sdg q16[1];
h q16[1];
h q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
h q16[0];
h q16[1];
s q16[1];
z q16[1];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
h q16[1];
s q16[1];
y q16[1];
s q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
y q16[0];
s q16[0];
x q16[1];
cx q16[0],q16[1];
y q16[0];
y q16[1];
s q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
x q16[0];
h q16[0];
h q16[1];
s q16[1];
z q16[1];
s q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
h q16[1];
s q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
z q16[0];
sdg q16[0];
h q16[0];
h q16[1];
cx q16[0],q16[1];
x q16[0];
sdg q16[0];
h q16[0];
h q16[1];
s q16[1];
x q16[1];
h q16[1];
cx q16[0],q16[1];
x q16[0];
h q16[0];
s q16[0];
sdg q16[1];
h q16[1];
y q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
y q16[0];
sdg q16[0];
h q16[0];
h q16[1];
s q16[1];
x q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
h q16[1];
s q16[1];
cx q16[0],q16[1];
x q16[0];
sdg q16[0];
h q16[0];
z q16[1];
s q16[1];
h q16[1];
s q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
h q16[1];
s q16[1];
cx q16[0],q16[1];
z q16[0];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
h q16[1];
s q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
x q16[0];
h q16[0];
s q16[0];
z q16[1];
s q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
h q16[1];
s q16[1];
cx q16[0],q16[1];
x q16[0];
s q16[0];
h q16[0];
s q16[0];
y q16[1];
s q16[1];
h q16[1];
s q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
z q16[0];
sdg q16[1];
h q16[1];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
z q16[0];
x q16[1];
s q16[1];
h q16[1];
s q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
y q16[0];
h q16[0];
y q16[1];
h q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
x q16[0];
s q16[0];
h q16[0];
s q16[0];
sdg q16[1];
h q16[1];
y q16[1];
h q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
h q16[1];
s q16[1];
cx q16[0],q16[1];
s q16[0];
z q16[1];
s q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
y q16[0];
h q16[0];
sdg q16[1];
h q16[1];
s q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
s q16[0];
h q16[0];
s q16[0];
z q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
z q16[0];
s q16[0];
h q16[1];
s q16[1];
y q16[1];
s q16[1];
h q16[1];
s q16[1];
cx q16[0],q16[1];
y q16[0];
s q16[0];
sdg q16[1];
h q16[1];
h q16[1];
s q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
x q16[0];
h q16[1];
s q16[1];
y q16[1];
s q16[1];
h q16[1];
s q16[1];
cx q16[0],q16[1];
x q16[0];
x q16[1];
h q16[1];
s q16[1];
cx q16[0],q16[1];
z q16[0];
sdg q16[0];
h q16[0];
x q16[1];
h q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
h q16[1];
s q16[1];
cx q16[0],q16[1];
z q16[0];
s q16[0];
h q16[0];
s q16[0];
x q16[0];
s q16[0];
z q16[1];
h q16[1];
s q16[1];
x q16[1];
s q16[1];
h q16[1];
s q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
x q16[0];
s q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
z q16[0];
s q16[0];
h q16[0];
s q16[0];
h q16[1];
s q16[1];
y q16[1];
h q16[1];
cx q16[0],q16[1];
y q16[0];
h q16[0];
y q16[1];
h q16[1];
s q16[1];
cx q16[0],q16[1];
h q16[0];
h q16[1];
s q16[1];
y q16[1];
s q16[1];
h q16[1];
s q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
h q16[0];
h q16[1];
s q16[1];
x q16[1];
s q16[1];
h q16[1];
s q16[1];
cx q16[0],q16[1];
z q16[0];
s q16[0];
h q16[0];
s q16[0];
sdg q16[1];
h q16[1];
z q16[1];
s q16[1];
h q16[1];
s q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
z q16[0];
h q16[0];
s q16[0];
h q16[1];
s q16[1];
y q16[1];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
h q16[1];
s q16[1];
cx q16[0],q16[1];
z q16[0];
s q16[0];
h q16[0];
s q16[0];
x q16[1];
h q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
x q16[0];
s q16[0];
h q16[1];
s q16[1];
x q16[1];
h q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
y q16[0];
s q16[0];
x q16[1];
h q16[1];
s q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
h q16[1];
s q16[1];
cx q16[0],q16[1];
z q16[0];
h q16[0];
s q16[0];
h q16[1];
s q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
z q16[0];
s q16[0];
sdg q16[1];
h q16[1];
h q16[1];
s q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
s q16[0];
h q16[0];
s q16[0];
x q16[1];
s q16[1];
h q16[1];
s q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
x q16[0];
h q16[0];
s q16[0];
y q16[1];
x q16[1];
s q16[1];
h q16[1];
s q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
x q16[0];
h q16[0];
s q16[0];
h q16[1];
s q16[1];
z q16[1];
h q16[1];
cx q16[0],q16[1];
y q16[0];
s q16[0];
h q16[0];
s q16[0];
x q16[0];
h q16[0];
s q16[0];
x q16[1];
z q16[1];
s q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
y q16[0];
y q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
z q16[0];
h q16[0];
s q16[0];
y q16[1];
s q16[1];
h q16[1];
s q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
z q16[0];
s q16[1];
cx q16[0],q16[1];
z q16[0];
s q16[0];
h q16[0];
s q16[0];
s q16[1];
cx q16[0],q16[1];
x q16[0];
h q16[0];
s q16[0];
h q16[1];
s q16[1];
h q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
h q16[0];
s q16[0];
h q16[1];
s q16[1];
h q16[1];
s q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
y q16[0];
s q16[0];
x q16[1];
s q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
h q16[0];
y q16[1];
s q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
y q16[0];
s q16[0];
sdg q16[1];
h q16[1];
y q16[1];
h q16[1];
s q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
y q16[0];
h q16[0];
h q16[1];
s q16[1];
z q16[1];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
x q16[0];
h q16[1];
s q16[1];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
h q16[1];
s q16[1];
cx q16[0],q16[1];
y q16[0];
s q16[0];
h q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
h q16[1];
s q16[1];
cx q16[0],q16[1];
z q16[0];
sdg q16[0];
h q16[0];
x q16[1];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
z q16[0];
z q16[1];
h q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
x q16[0];
sdg q16[0];
h q16[0];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
y q16[0];
h q16[0];
sdg q16[1];
h q16[1];
x q16[1];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
s q16[0];
h q16[0];
s q16[0];
sdg q16[1];
h q16[1];
z q16[1];
h q16[1];
s q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
h q16[1];
s q16[1];
cx q16[0],q16[1];
y q16[0];
s q16[0];
z q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
z q16[0];
sdg q16[0];
h q16[0];
h q16[1];
s q16[1];
y q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
h q16[0];
z q16[1];
h q16[1];
cx q16[0],q16[1];
z q16[0];
h q16[1];
s q16[1];
y q16[1];
h q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
y q16[0];
sdg q16[1];
h q16[1];
x q16[1];
s q16[1];
h q16[1];
s q16[1];
cx q16[0],q16[1];
s q16[0];
h q16[0];
s q16[0];
x q16[1];
h q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
h q16[1];
s q16[1];
z q16[1];
s q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
s q16[0];
h q16[1];
s q16[1];
y q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
y q16[0];
sdg q16[0];
h q16[0];
h q16[1];
s q16[1];
h q16[1];
s q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
y q16[0];
sdg q16[0];
h q16[0];
s q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
h q16[1];
s q16[1];
cx q16[0],q16[1];
z q16[0];
s q16[0];
h q16[0];
s q16[0];
y q16[1];
s q16[1];
cx q16[0],q16[1];
y q16[0];
sdg q16[1];
h q16[1];
z q16[1];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
z q16[0];
s q16[0];
h q16[0];
s q16[0];
sdg q16[1];
h q16[1];
y q16[1];
h q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
x q16[0];
s q16[0];
x q16[1];
h q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
h q16[1];
s q16[1];
cx q16[0],q16[1];
s q16[0];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
z q16[1];
s q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
h q16[0];
z q16[1];
s q16[1];
h q16[1];
s q16[1];
cx q16[0],q16[1];
y q16[0];
s q16[0];
h q16[0];
s q16[0];
z q16[1];
cx q16[0],q16[1];
x q16[0];
s q16[0];
h q16[0];
s q16[0];
z q16[1];
s q16[1];
h q16[1];
s q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
x q16[0];
h q16[0];
s q16[0];
h q16[1];
s q16[1];
s q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
h q16[1];
s q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
z q16[0];
h q16[0];
x q16[1];
s q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
z q16[0];
h q16[0];
s q16[0];
h q16[1];
s q16[1];
y q16[1];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
z q16[0];
h q16[0];
s q16[0];
h q16[1];
s q16[1];
z q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
h q16[1];
s q16[1];
cx q16[0],q16[1];
y q16[0];
h q16[0];
x q16[1];
s q16[1];
h q16[1];
s q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
h q16[1];
s q16[1];
cx q16[0],q16[1];
z q16[0];
s q16[0];
h q16[0];
s q16[0];
z q16[1];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
h q16[1];
s q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
s q16[0];
h q16[0];
s q16[0];
z q16[1];
s q16[1];
h q16[1];
s q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
y q16[0];
x q16[1];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
z q16[0];
h q16[0];
s q16[0];
x q16[1];
h q16[1];
cx q16[0],q16[1];
y q16[0];
s q16[0];
h q16[0];
s q16[0];
sdg q16[1];
h q16[1];
z q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
h q16[1];
s q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
x q16[0];
s q16[0];
x q16[1];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
s q16[0];
h q16[0];
s q16[0];
z q16[1];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
x q16[0];
h q16[0];
h q16[1];
s q16[1];
s q16[1];
h q16[1];
s q16[1];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
sdg q16[1];
h q16[1];
z q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
x q16[0];
s q16[0];
h q16[0];
s q16[0];
h q16[1];
s q16[1];
z q16[1];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
h q16[1];
s q16[1];
cx q16[0],q16[1];
z q16[0];
x q16[1];
s q16[1];
h q16[1];
s q16[1];
cx q16[0],q16[1];
y q16[0];
h q16[0];
h q16[1];
s q16[1];
z q16[1];
s q16[1];
h q16[1];
s q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
z q16[0];
sdg q16[0];
h q16[0];
h q16[1];
s q16[1];
x q16[1];
h q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
y q16[0];
sdg q16[1];
h q16[1];
y q16[1];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
y q16[0];
y q16[1];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
s q16[0];
h q16[1];
s q16[1];
y q16[1];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
x q16[0];
h q16[0];
sdg q16[1];
h q16[1];
h q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
x q16[0];
h q16[0];
x q16[1];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
h q16[0];
s q16[0];
z q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
x q16[0];
s q16[0];
cx q16[0],q16[1];
h q16[0];
s q16[0];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
x q16[0];
s q16[0];
x q16[1];
s q16[1];
h q16[1];
s q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
h q16[1];
s q16[1];
x q16[1];
s q16[1];
h q16[1];
s q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
z q16[1];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
s q16[0];
h q16[0];
s q16[0];
s q16[0];
x q16[1];
sdg q16[1];
h q16[1];
z q16[1];
h q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
h q16[1];
s q16[1];
cx q16[0],q16[1];
x q16[0];
h q16[0];
cx q16[0],q16[1];
h q16[0];
s q16[0];
h q16[1];
s q16[1];
cx q16[0],q16[1];
z q16[0];
y q16[1];
cx q16[0],q16[1];
x q16[0];
h q16[0];
y q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
y q16[0];
s q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
s q16[0];
h q16[0];
s q16[0];
sdg q16[1];
h q16[1];
x q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
x q16[0];
z q16[1];
cx q16[0],q16[1];
y q16[0];
s q16[0];
h q16[1];
s q16[1];
y q16[1];
s q16[1];
h q16[1];
s q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
y q16[0];
z q16[1];
s q16[1];
h q16[1];
s q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
x q16[0];
h q16[0];
y q16[0];
s q16[0];
x q16[0];
s q16[0];
h q16[0];
s q16[0];
s q16[1];
y q16[1];
x q16[1];
s q16[1];
h q16[1];
s q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
y q16[0];
sdg q16[0];
h q16[0];
h q16[1];
s q16[1];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
x q16[0];
s q16[0];
x q16[0];
h q16[1];
s q16[1];
y q16[1];
s q16[1];
h q16[1];
s q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
y q16[1];
h q16[1];
s q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
y q16[0];
h q16[0];
s q16[0];
s q16[1];
h q16[1];
s q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
h q16[1];
s q16[1];
cx q16[0],q16[1];
y q16[1];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
x q16[0];
sdg q16[0];
h q16[0];
s q16[1];
h q16[1];
s q16[1];
cx q16[0],q16[1];
z q16[0];
sdg q16[0];
h q16[0];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
h q16[1];
s q16[1];
cx q16[0],q16[1];
z q16[0];
sdg q16[0];
h q16[0];
y q16[1];
h q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
z q16[0];
s q16[0];
h q16[0];
s q16[0];
h q16[1];
s q16[1];
z q16[1];
cx q16[0],q16[1];
z q16[0];
h q16[0];
s q16[0];
sdg q16[1];
h q16[1];
h q16[1];
x q16[1];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
h q16[1];
s q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
y q16[0];
x q16[1];
s q16[1];
h q16[1];
s q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
z q16[0];
sdg q16[0];
h q16[0];
y q16[1];
s q16[1];
h q16[1];
s q16[1];
cx q16[0],q16[1];
x q16[0];
sdg q16[0];
h q16[0];
sdg q16[1];
h q16[1];
x q16[1];
h q16[1];
s q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
z q16[0];
h q16[0];
h q16[1];
s q16[1];
x q16[1];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
y q16[0];
sdg q16[0];
h q16[0];
h q16[1];
s q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
y q16[0];
h q16[1];
s q16[1];
z q16[1];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
z q16[0];
h q16[1];
s q16[1];
z q16[1];
s q16[1];
h q16[1];
s q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
h q16[1];
s q16[1];
cx q16[0],q16[1];
s q16[0];
s q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
y q16[1];
h q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
sdg q16[1];
h q16[1];
x q16[1];
s q16[1];
h q16[1];
s q16[1];
cx q16[0],q16[1];
s q16[1];
h q16[1];
s q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
z q16[0];
s q16[0];
h q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
y q16[0];
h q16[0];
s q16[0];
h q16[1];
s q16[1];
y q16[1];
h q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
y q16[0];
s q16[0];
y q16[1];
s q16[1];
h q16[1];
s q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
x q16[0];
h q16[0];
s q16[0];
h q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
h q16[1];
s q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
s q16[0];
h q16[0];
s q16[0];
x q16[1];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
y q16[0];
sdg q16[0];
h q16[0];
z q16[1];
s q16[1];
h q16[1];
s q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
y q16[0];
s q16[0];
h q16[0];
s q16[0];
z q16[1];
s q16[1];
h q16[1];
s q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
y q16[0];
s q16[0];
h q16[0];
s q16[0];
x q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
h q16[1];
s q16[1];
cx q16[0],q16[1];
y q16[0];
h q16[0];
s q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
h q16[1];
s q16[1];
cx q16[0],q16[1];
h q16[0];
z q16[1];
h q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
x q16[0];
s q16[0];
h q16[1];
s q16[1];
z q16[1];
s q16[1];
h q16[1];
s q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
x q16[0];
sdg q16[0];
h q16[0];
y q16[1];
s q16[1];
h q16[1];
s q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
x q16[0];
s q16[0];
h q16[0];
s q16[0];
h q16[1];
s q16[1];
z q16[1];
s q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
h q16[1];
s q16[1];
cx q16[0],q16[1];
x q16[1];
h q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
s q16[0];
y q16[1];
s q16[1];
cx q16[0],q16[1];
y q16[0];
s q16[0];
h q16[1];
s q16[1];
x q16[1];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
x q16[0];
h q16[1];
s q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
x q16[0];
h q16[0];
h q16[1];
s q16[1];
x q16[1];
s q16[1];
h q16[1];
s q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
x q16[0];
s q16[0];
h q16[0];
s q16[0];
x q16[1];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
x q16[0];
s q16[0];
h q16[1];
s q16[1];
x q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
x q16[0];
s q16[0];
h q16[1];
s q16[1];
x q16[1];
s q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
x q16[0];
x q16[1];
s q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
x q16[0];
h q16[1];
s q16[1];
x q16[1];
h q16[1];
s q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
h q16[1];
s q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
x q16[0];
h q16[0];
s q16[0];
y q16[1];
s q16[1];
h q16[1];
s q16[1];
cx q16[0],q16[1];
s q16[0];
x q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
s q16[0];
h q16[0];
s q16[0];
h q16[1];
s q16[1];
z q16[1];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
z q16[0];
sdg q16[0];
h q16[0];
x q16[1];
h q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
h q16[1];
s q16[1];
cx q16[0],q16[1];
y q16[0];
h q16[0];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
h q16[1];
s q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
x q16[0];
sdg q16[0];
h q16[0];
z q16[1];
s q16[1];
h q16[1];
s q16[1];
x q16[1];
s q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
h q16[1];
s q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
y q16[0];
s q16[0];
h q16[0];
s q16[0];
h q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
x q16[0];
y q16[1];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
y q16[0];
h q16[0];
h q16[1];
s q16[1];
h q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
y q16[0];
sdg q16[0];
h q16[0];
z q16[1];
s q16[1];
h q16[1];
s q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
h q16[0];
s q16[0];
z q16[0];
sdg q16[0];
h q16[0];
y q16[0];
h q16[0];
s q16[0];
z q16[1];
h q16[1];
x q16[1];
s q16[1];
h q16[1];
s q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
y q16[0];
h q16[0];
h q16[1];
s q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
z q16[1];
s q16[1];
h q16[1];
s q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
h q16[1];
s q16[1];
cx q16[0],q16[1];
y q16[0];
h q16[0];
s q16[0];
cx q16[0],q16[1];
x q16[0];
sdg q16[1];
h q16[1];
y q16[1];
s q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
y q16[0];
s q16[0];
h q16[0];
s q16[0];
h q16[1];
s q16[1];
z q16[1];
h q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
s q16[0];
h q16[0];
s q16[0];
y q16[1];
h q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
h q16[1];
s q16[1];
cx q16[0],q16[1];
y q16[0];
h q16[0];
z q16[1];
s q16[1];
h q16[1];
s q16[1];
cx q16[0],q16[1];
y q16[0];
h q16[0];
s q16[0];
h q16[1];
s q16[1];
z q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
y q16[0];
h q16[0];
h q16[1];
s q16[1];
h q16[1];
s q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
h q16[1];
s q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
x q16[1];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
z q16[0];
s q16[0];
sdg q16[1];
h q16[1];
z q16[1];
s q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
sdg q16[0];
h q16[0];
sdg q16[1];
h q16[1];
x q16[1];
h q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
h q16[1];
s q16[1];
cx q16[0],q16[1];
s q16[0];
h q16[0];
s q16[0];
z q16[0];
h q16[0];
s q16[0];
z q16[1];
s q16[1];
h q16[1];
s q16[1];
z q16[1];
h q16[1];
s q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
h q16[1];
s q16[1];
cx q16[0],q16[1];
y q16[0];
z q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
x q16[0];
h q16[0];
h q16[1];
s q16[1];
y q16[1];
h q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
h q16[1];
s q16[1];
cx q16[0],q16[1];
x q16[0];
s q16[0];
h q16[0];
s q16[0];
y q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
h q16[1];
s q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
x q16[0];
sdg q16[0];
h q16[0];
y q16[1];
h q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
y q16[0];
h q16[0];
s q16[0];
z q16[1];
s q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
h q16[1];
s q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
z q16[0];
s q16[0];
h q16[0];
s q16[0];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
y q16[0];
sdg q16[0];
h q16[0];
y q16[1];
s q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
h q16[1];
s q16[1];
cx q16[0],q16[1];
x q16[0];
s q16[0];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
z q16[0];
s q16[0];
h q16[0];
s q16[0];
x q16[1];
h q16[1];
s q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
h q16[1];
s q16[1];
cx q16[0],q16[1];
z q16[0];
s q16[0];
x q16[1];
h q16[1];
s q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
z q16[0];
h q16[0];
s q16[0];
h q16[1];
s q16[1];
y q16[1];
s q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
h q16[1];
s q16[1];
s q16[1];
h q16[1];
s q16[1];
cx q16[0],q16[1];
h q16[1];
s q16[1];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
h q16[1];
s q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
x q16[1];
s q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
x q16[0];
h q16[0];
s q16[0];
sdg q16[1];
h q16[1];
x q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
y q16[0];
s q16[0];
h q16[1];
s q16[1];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
x q16[0];
s q16[0];
sdg q16[1];
h q16[1];
z q16[1];
h q16[1];
cx q16[0],q16[1];
z q16[0];
h q16[0];
s q16[0];
h q16[1];
s q16[1];
y q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
x q16[0];
s q16[0];
h q16[0];
s q16[0];
h q16[1];
s q16[1];
y q16[1];
s q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
h q16[1];
s q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
z q16[0];
h q16[0];
y q16[1];
h q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
s q16[0];
h q16[0];
s q16[0];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
y q16[0];
h q16[0];
s q16[0];
sdg q16[1];
h q16[1];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
h q16[1];
s q16[1];
cx q16[0],q16[1];
x q16[0];
s q16[0];
y q16[0];
s q16[0];
h q16[0];
s q16[0];
x q16[1];
h q16[1];
y q16[1];
s q16[1];
h q16[1];
s q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
h q16[1];
s q16[1];
cx q16[0],q16[1];
x q16[0];
h q16[0];
s q16[0];
y q16[1];
h q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
h q16[1];
s q16[1];
cx q16[0],q16[1];
x q16[0];
x q16[1];
s q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
x q16[0];
s q16[0];
s q16[0];
h q16[0];
s q16[0];
z q16[1];
h q16[1];
s q16[1];
y q16[1];
s q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
h q16[1];
s q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
z q16[0];
h q16[0];
s q16[0];
x q16[1];
s q16[1];
h q16[1];
s q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
h q16[1];
s q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
z q16[1];
s q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
y q16[0];
s q16[0];
x q16[1];
h q16[1];
s q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
y q16[0];
x q16[1];
h q16[1];
s q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
s q16[0];
h q16[0];
s q16[0];
h q16[1];
s q16[1];
y q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
x q16[0];
s q16[0];
h q16[0];
s q16[0];
y q16[1];
s q16[1];
h q16[1];
s q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
h q16[1];
s q16[1];
cx q16[0],q16[1];
y q16[0];
h q16[0];
y q16[1];
h q16[1];
s q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
h q16[0];
x q16[1];
h q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
h q16[1];
s q16[1];
cx q16[0],q16[1];
z q16[0];
h q16[0];
s q16[0];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
y q16[0];
h q16[1];
s q16[1];
s q16[1];
h q16[1];
s q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
h q16[1];
s q16[1];
cx q16[0],q16[1];
z q16[0];
s q16[0];
y q16[0];
h q16[0];
s q16[0];
sdg q16[1];
h q16[1];
z q16[1];
h q16[1];
cx q16[0],q16[1];
z q16[0];
sdg q16[1];
h q16[1];
x q16[1];
s q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
y q16[1];
s q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
z q16[0];
x q16[1];
s q16[1];
cx q16[0],q16[1];
z q16[0];
sdg q16[0];
h q16[0];
h q16[1];
s q16[1];
x q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
x q16[0];
s q16[0];
sdg q16[1];
h q16[1];
y q16[1];
h q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
z q16[0];
s q16[0];
h q16[0];
s q16[0];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
y q16[0];
s q16[0];
x q16[1];
h q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
z q16[0];
sdg q16[0];
h q16[0];
sdg q16[1];
h q16[1];
x q16[1];
s q16[1];
h q16[1];
s q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
z q16[0];
s q16[0];
h q16[1];
s q16[1];
z q16[1];
h q16[1];
s q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
sdg q16[0];
h q16[0];
s q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
h q16[1];
s q16[1];
x q16[1];
s q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
z q16[0];
s q16[0];
y q16[1];
s q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
y q16[0];
h q16[1];
s q16[1];
h q16[1];
s q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
z q16[0];
s q16[0];
h q16[0];
s q16[0];
h q16[1];
s q16[1];
y q16[1];
h q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
y q16[1];
h q16[1];
s q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
x q16[0];
h q16[0];
h q16[1];
s q16[1];
x q16[1];
h q16[1];
s q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
h q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
y q16[0];
h q16[0];
y q16[1];
h q16[1];
s q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
z q16[0];
sdg q16[0];
h q16[0];
h q16[1];
s q16[1];
z q16[1];
s q16[1];
h q16[1];
s q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
x q16[0];
s q16[0];
h q16[0];
s q16[0];
h q16[1];
s q16[1];
y q16[1];
h q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
z q16[0];
sdg q16[0];
h q16[0];
h q16[1];
s q16[1];
z q16[1];
s q16[1];
h q16[1];
s q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
h q16[1];
s q16[1];
h q16[1];
s q16[1];
cx q16[0],q16[1];
y q16[0];
s q16[0];
h q16[1];
s q16[1];
s q16[1];
h q16[1];
s q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
z q16[0];
h q16[0];
s q16[0];
x q16[1];
cx q16[0],q16[1];
x q16[0];
sdg q16[1];
h q16[1];
z q16[1];
s q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
z q16[1];
s q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
h q16[1];
s q16[1];
cx q16[0],q16[1];
h q16[0];
y q16[1];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
y q16[0];
s q16[0];
s q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
z q16[0];
sdg q16[0];
h q16[0];
z q16[1];
s q16[1];
h q16[1];
s q16[1];
cx q16[0],q16[1];
y q16[0];
h q16[0];
s q16[0];
sdg q16[1];
h q16[1];
y q16[1];
h q16[1];
s q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
x q16[0];
h q16[0];
s q16[0];
h q16[1];
s q16[1];
y q16[1];
s q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
s q16[0];
h q16[0];
s q16[0];
sdg q16[1];
h q16[1];
z q16[1];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
h q16[1];
s q16[1];
y q16[1];
h q16[1];
s q16[1];
cx q16[0],q16[1];
y q16[0];
h q16[0];
s q16[0];
h q16[1];
s q16[1];
cx q16[0],q16[1];
s q16[0];
h q16[0];
s q16[0];
h q16[1];
s q16[1];
z q16[1];
h q16[1];
cx q16[0],q16[1];
s q16[0];
sdg q16[1];
h q16[1];
s q16[1];
cx q16[0],q16[1];
y q16[0];
s q16[0];
sdg q16[1];
h q16[1];
x q16[1];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
h q16[0];
s q16[0];
sdg q16[1];
h q16[1];
cx q16[0],q16[1];
x q16[0];
s q16[0];
h q16[0];
s q16[0];
h q16[1];
s q16[1];
z q16[1];
s q16[1];
h q16[1];
s q16[1];
cx q16[0],q16[1];
sdg q16[0];
h q16[0];
h q16[1];
s q16[1];
