OPENQASM 2.0;
include "qelib1.inc";
qreg q38[2];
x q38[0];
s q38[0];
z q38[1];
cx q38[0],q38[1];
sdg q38[0];
h q38[0];
sdg q38[1];
h q38[1];
cx q38[0],q38[1];
h q38[0];
s q38[0];
z q38[0];
s q38[0];
h q38[0];
s q38[0];
h q38[1];
s q38[1];
cx q38[0],q38[1];
h q38[0];
s q38[0];
sdg q38[1];
h q38[1];
cx q38[0],q38[1];
z q38[0];
h q38[1];
s q38[1];
x q38[1];
s q38[1];
cx q38[0],q38[1];
sdg q38[0];
h q38[0];
sdg q38[1];
h q38[1];
cx q38[0],q38[1];
x q38[0];
h q38[0];
s q38[1];
h q38[1];
s q38[1];
cx q38[0],q38[1];
sdg q38[0];
h q38[0];
x q38[0];
h q38[0];
s q38[0];
sdg q38[1];
h q38[1];
z q38[1];
cx q38[0],q38[1];
sdg q38[0];
h q38[0];
x q38[1];
sdg q38[1];
h q38[1];
cx q38[0],q38[1];
sdg q38[0];
h q38[0];
sdg q38[1];
h q38[1];
cx q38[0],q38[1];
h q38[0];
s q38[0];
h q38[1];
s q38[1];
cx q38[0],q38[1];
s q38[0];
cx q38[0],q38[1];
h q38[0];
s q38[0];
z q38[0];
sdg q38[0];
h q38[0];
h q38[1];
s q38[1];
z q38[1];
sdg q38[1];
h q38[1];
cx q38[0],q38[1];
sdg q38[0];
h q38[0];
sdg q38[1];
h q38[1];
cx q38[0],q38[1];
x q38[0];
sdg q38[0];
h q38[0];
x q38[1];
s q38[1];
h q38[1];
s q38[1];
cx q38[0],q38[1];
sdg q38[0];
h q38[0];
sdg q38[1];
h q38[1];
cx q38[0],q38[1];
s q38[0];
h q38[0];
s q38[0];
z q38[1];
sdg q38[1];
h q38[1];
cx q38[0],q38[1];
sdg q38[0];
h q38[0];
z q38[0];
s q38[0];
h q38[0];
s q38[0];
x q38[1];
cx q38[0],q38[1];
sdg q38[0];
h q38[0];
x q38[0];
sdg q38[1];
h q38[1];
cx q38[0],q38[1];
h q38[0];
s q38[0];
x q38[0];
s q38[0];
h q38[1];
s q38[1];
z q38[1];
h q38[1];
s q38[1];
cx q38[0],q38[1];
x q38[0];
s q38[0];
h q38[0];
s q38[0];
z q38[1];
s q38[1];
cx q38[0],q38[1];
sdg q38[0];
h q38[0];
sdg q38[1];
h q38[1];
cx q38[0],q38[1];
x q38[0];
s q38[0];
h q38[1];
s q38[1];
sdg q38[1];
h q38[1];
cx q38[0],q38[1];
sdg q38[0];
h q38[0];
sdg q38[1];
h q38[1];
cx q38[0],q38[1];
y q38[0];
h q38[0];
h q38[1];
s q38[1];
cx q38[0],q38[1];
x q38[0];
h q38[0];
s q38[0];
h q38[1];
s q38[1];
x q38[1];
sdg q38[1];
h q38[1];
cx q38[0],q38[1];
h q38[0];
s q38[0];
sdg q38[1];
h q38[1];
cx q38[0],q38[1];
y q38[0];
s q38[0];
h q38[0];
s q38[0];
h q38[1];
s q38[1];
y q38[1];
h q38[1];
s q38[1];
cx q38[0],q38[1];
sdg q38[1];
h q38[1];
z q38[1];
s q38[1];
cx q38[0],q38[1];
sdg q38[0];
h q38[0];
h q38[1];
s q38[1];
cx q38[0],q38[1];
x q38[0];
s q38[0];
z q38[1];
cx q38[0],q38[1];
sdg q38[0];
h q38[0];
sdg q38[1];
h q38[1];
cx q38[0],q38[1];
h q38[0];
s q38[0];
z q38[0];
s q38[0];
z q38[1];
h q38[1];
s q38[1];
cx q38[0],q38[1];
h q38[0];
s q38[0];
y q38[0];
s q38[0];
h q38[1];
s q38[1];
x q38[1];
h q38[1];
s q38[1];
cx q38[0],q38[1];
h q38[0];
s q38[0];
y q38[0];
sdg q38[0];
h q38[0];
sdg q38[1];
h q38[1];
cx q38[0],q38[1];
x q38[0];
h q38[0];
sdg q38[1];
h q38[1];
x q38[1];
s q38[1];
h q38[1];
s q38[1];
cx q38[0],q38[1];
sdg q38[0];
h q38[0];
h q38[1];
s q38[1];
cx q38[0],q38[1];
h q38[0];
s q38[0];
x q38[0];
s q38[0];
z q38[1];
h q38[1];
s q38[1];
cx q38[0],q38[1];
z q38[0];
h q38[0];
sdg q38[1];
h q38[1];
z q38[1];
h q38[1];
cx q38[0],q38[1];
sdg q38[0];
h q38[0];
sdg q38[1];
h q38[1];
cx q38[0],q38[1];
z q38[0];
s q38[0];
h q38[0];
s q38[0];
h q38[1];
s q38[1];
y q38[1];
h q38[1];
s q38[1];
cx q38[0],q38[1];
h q38[0];
s q38[0];
sdg q38[1];
h q38[1];
cx q38[0],q38[1];
y q38[0];
h q38[1];
s q38[1];
x q38[1];
h q38[1];
s q38[1];
cx q38[0],q38[1];
z q38[0];
sdg q38[0];
h q38[0];
x q38[1];
cx q38[0],q38[1];
h q38[0];
s q38[0];
y q38[0];
sdg q38[1];
h q38[1];
y q38[1];
h q38[1];
s q38[1];
cx q38[0],q38[1];
h q38[0];
s q38[0];
h q38[1];
s q38[1];
cx q38[0],q38[1];
y q38[0];
h q38[0];
s q38[0];
h q38[1];
cx q38[0],q38[1];
sdg q38[0];
h q38[0];
sdg q38[1];
h q38[1];
cx q38[0],q38[1];
s q38[0];
h q38[0];
s q38[0];
h q38[1];
s q38[1];
x q38[1];
cx q38[0],q38[1];
sdg q38[0];
h q38[0];
z q38[0];
h q38[0];
sdg q38[1];
h q38[1];
cx q38[0],q38[1];
sdg q38[0];
h q38[0];
h q38[1];
s q38[1];
cx q38[0],q38[1];
h q38[0];
s q38[0];
sdg q38[0];
h q38[0];
z q38[1];
s q38[1];
cx q38[0],q38[1];
sdg q38[0];
h q38[0];
sdg q38[1];
h q38[1];
cx q38[0],q38[1];
h q38[0];
s q38[0];
h q38[1];
s q38[1];
cx q38[0],q38[1];
x q38[0];
s q38[0];
x q38[1];
s q38[1];
h q38[1];
s q38[1];
cx q38[0],q38[1];
h q38[0];
s q38[0];
sdg q38[1];
h q38[1];
cx q38[0],q38[1];
z q38[0];
h q38[0];
s q38[0];
h q38[1];
s q38[1];
z q38[1];
cx q38[0],q38[1];
sdg q38[0];
h q38[0];
y q38[0];
h q38[0];
h q38[1];
s q38[1];
y q38[1];
cx q38[0],q38[1];
h q38[0];
s q38[0];
sdg q38[1];
h q38[1];
cx q38[0],q38[1];
z q38[0];
h q38[0];
h q38[1];
s q38[1];
x q38[1];
cx q38[0],q38[1];
sdg q38[0];
h q38[0];
h q38[1];
s q38[1];
cx q38[0],q38[1];
z q38[0];
s q38[0];
y q38[1];
sdg q38[1];
h q38[1];
cx q38[0],q38[1];
sdg q38[0];
h q38[0];
h q38[0];
h q38[1];
s q38[1];
y q38[1];
s q38[1];
cx q38[0],q38[1];
z q38[0];
sdg q38[0];
h q38[0];
sdg q38[1];
h q38[1];
x q38[1];
s q38[1];
h q38[1];
s q38[1];
cx q38[0],q38[1];
sdg q38[0];
h q38[0];
h q38[1];
s q38[1];
cx q38[0],q38[1];
h q38[0];
s q38[0];
y q38[0];
h q38[0];
y q38[1];
cx q38[0],q38[1];
sdg q38[0];
h q38[0];
sdg q38[1];
h q38[1];
cx q38[0],q38[1];
s q38[0];
h q38[0];
s q38[0];
y q38[1];
h q38[1];
s q38[1];
cx q38[0],q38[1];
sdg q38[0];
h q38[0];
sdg q38[1];
h q38[1];
cx q38[0],q38[1];
h q38[0];
s q38[0];
z q38[0];
h q38[0];
h q38[1];
s q38[1];
z q38[1];
sdg q38[1];
h q38[1];
cx q38[0],q38[1];
sdg q38[0];
h q38[0];
sdg q38[1];
h q38[1];
cx q38[0],q38[1];
z q38[0];
h q38[0];
y q38[1];
s q38[1];
h q38[1];
s q38[1];
cx q38[0],q38[1];
h q38[0];
s q38[0];
h q38[1];
s q38[1];
cx q38[0],q38[1];
x q38[0];
s q38[0];
h q38[0];
s q38[0];
x q38[1];
cx q38[0],q38[1];
h q38[0];
s q38[0];
z q38[0];
sdg q38[0];
h q38[0];
h q38[1];
s q38[1];
x q38[1];
sdg q38[1];
h q38[1];
cx q38[0],q38[1];
h q38[0];
s q38[0];
h q38[1];
s q38[1];
cx q38[0],q38[1];
z q38[0];
sdg q38[0];
h q38[0];
y q38[1];
h q38[1];
cx q38[0],q38[1];
h q38[0];
s q38[0];
h q38[1];
s q38[1];
cx q38[0],q38[1];
s q38[0];
h q38[0];
s q38[0];
x q38[1];
h q38[1];
cx q38[0],q38[1];
sdg q38[0];
h q38[0];
h q38[1];
s q38[1];
cx q38[0],q38[1];
h q38[0];
s q38[0];
x q38[1];
sdg q38[1];
h q38[1];
cx q38[0],q38[1];
h q38[0];
s q38[0];
y q38[0];
sdg q38[1];
h q38[1];
h q38[1];
cx q38[0],q38[1];
h q38[0];
s q38[0];
h q38[1];
s q38[1];
cx q38[0],q38[1];
h q38[0];
s q38[0];
x q38[1];
h q38[1];
s q38[1];
cx q38[0],q38[1];
s q38[0];
h q38[0];
s q38[0];
h q38[1];
s q38[1];
z q38[1];
h q38[1];
s q38[1];
cx q38[0],q38[1];
sdg q38[0];
h q38[0];
h q38[1];
s q38[1];
cx q38[0],q38[1];
h q38[0];
s q38[0];
y q38[0];
h q38[0];
s q38[0];
h q38[1];
cx q38[0],q38[1];
sdg q38[0];
h q38[0];
sdg q38[1];
h q38[1];
cx q38[0],q38[1];
h q38[0];
s q38[0];
z q38[0];
h q38[0];
s q38[0];
h q38[1];
s q38[1];
z q38[1];
h q38[1];
s q38[1];
cx q38[0],q38[1];
x q38[0];
h q38[0];
s q38[0];
sdg q38[1];
h q38[1];
z q38[1];
h q38[1];
s q38[1];
cx q38[0],q38[1];
sdg q38[0];
h q38[0];
sdg q38[1];
h q38[1];
cx q38[0],q38[1];
z q38[0];
h q38[0];
s q38[0];
z q38[1];
h q38[1];
cx q38[0],q38[1];
sdg q38[0];
h q38[0];
sdg q38[1];
h q38[1];
cx q38[0],q38[1];
h q38[0];
s q38[0];
z q38[0];
s q38[0];
h q38[0];
s q38[0];
h q38[1];
s q38[1];
h q38[1];
s q38[1];
cx q38[0],q38[1];
sdg q38[0];
h q38[0];
y q38[0];
h q38[0];
s q38[0];
sdg q38[1];
h q38[1];
h q38[1];
s q38[1];
cx q38[0],q38[1];
sdg q38[0];
h q38[0];
h q38[1];
s q38[1];
cx q38[0],q38[1];
x q38[0];
h q38[0];
y q38[1];
h q38[1];
cx q38[0],q38[1];
h q38[0];
s q38[0];
x q38[0];
s q38[0];
h q38[0];
s q38[0];
h q38[1];
s q38[1];
y q38[1];
h q38[1];
s q38[1];
cx q38[0],q38[1];
sdg q38[0];
h q38[0];
sdg q38[1];
h q38[1];
cx q38[0],q38[1];
h q38[0];
s q38[0];
z q38[0];
s q38[0];
h q38[0];
s q38[0];
z q38[1];
sdg q38[1];
h q38[1];
cx q38[0],q38[1];
h q38[0];
s q38[0];
y q38[0];
h q38[1];
s q38[1];
x q38[1];
cx q38[0],q38[1];
x q38[0];
h q38[0];
s q38[0];
h q38[1];
s q38[1];
x q38[1];
s q38[1];
cx q38[0],q38[1];
sdg q38[0];
h q38[0];
h q38[1];
s q38[1];
cx q38[0],q38[1];
s q38[0];
z q38[1];
h q38[1];
cx q38[0],q38[1];
sdg q38[0];
h q38[0];
sdg q38[1];
h q38[1];
cx q38[0],q38[1];
h q38[0];
s q38[0];
h q38[1];
s q38[1];
cx q38[0],q38[1];
y q38[0];
h q38[0];
x q38[1];
s q38[1];
h q38[1];
s q38[1];
cx q38[0],q38[1];
sdg q38[0];
h q38[0];
h q38[1];
s q38[1];
cx q38[0],q38[1];
h q38[0];
s q38[0];
x q38[1];
h q38[1];
cx q38[0],q38[1];
h q38[0];
s q38[0];
y q38[0];
s q38[0];
h q38[1];
s q38[1];
x q38[1];
s q38[1];
cx q38[0],q38[1];
h q38[0];
s q38[0];
h q38[1];
s q38[1];
cx q38[0],q38[1];
x q38[0];
s q38[0];
h q38[0];
s q38[0];
x q38[1];
s q38[1];
cx q38[0],q38[1];
sdg q38[0];
h q38[0];
sdg q38[1];
h q38[1];
cx q38[0],q38[1];
s q38[0];
x q38[1];
s q38[1];
h q38[1];
s q38[1];
cx q38[0],q38[1];
h q38[0];
s q38[0];
x q38[0];
h q38[0];
s q38[0];
sdg q38[1];
h q38[1];
cx q38[0],q38[1];
z q38[0];
h q38[0];
sdg q38[1];
h q38[1];
x q38[1];
sdg q38[1];
h q38[1];
cx q38[0],q38[1];
sdg q38[0];
h q38[0];
sdg q38[1];
h q38[1];
cx q38[0],q38[1];
h q38[0];
s q38[0];
h q38[1];
s q38[1];
cx q38[0],q38[1];
x q38[0];
h q38[0];
y q38[1];
s q38[1];
h q38[1];
s q38[1];
cx q38[0],q38[1];
h q38[0];
s q38[0];
sdg q38[1];
h q38[1];
cx q38[0],q38[1];
sdg q38[0];
h q38[0];
h q38[1];
cx q38[0],q38[1];
sdg q38[0];
h q38[0];
sdg q38[1];
h q38[1];
cx q38[0],q38[1];
h q38[0];
s q38[0];
h q38[1];
s q38[1];
cx q38[0],q38[1];
y q38[0];
sdg q38[0];
h q38[0];
s q38[0];
h q38[0];
s q38[0];
h q38[1];
s q38[1];
s q38[1];
cx q38[0],q38[1];
sdg q38[0];
h q38[0];
sdg q38[1];
h q38[1];
cx q38[0],q38[1];
y q38[0];
s q38[0];
h q38[1];
s q38[1];
x q38[1];
s q38[1];
h q38[1];
s q38[1];
cx q38[0],q38[1];
h q38[0];
s q38[0];
sdg q38[1];
h q38[1];
cx q38[0],q38[1];
s q38[0];
h q38[1];
s q38[1];
z q38[1];
h q38[1];
s q38[1];
cx q38[0],q38[1];
h q38[0];
s q38[0];
h q38[1];
s q38[1];
cx q38[0],q38[1];
h q38[0];
y q38[1];
h q38[1];
s q38[1];
cx q38[0],q38[1];
sdg q38[0];
h q38[0];
sdg q38[1];
h q38[1];
cx q38[0],q38[1];
h q38[0];
s q38[0];
h q38[1];
s q38[1];
cx q38[0],q38[1];
z q38[0];
sdg q38[0];
h q38[0];
z q38[1];
s q38[1];
cx q38[0],q38[1];
h q38[0];
s q38[0];
sdg q38[1];
h q38[1];
cx q38[0],q38[1];
z q38[0];
h q38[0];
s q38[0];
h q38[1];
s q38[1];
y q38[1];
h q38[1];
cx q38[0],q38[1];
sdg q38[0];
h q38[0];
sdg q38[1];
h q38[1];
cx q38[0],q38[1];
y q38[0];
sdg q38[0];
h q38[0];
h q38[1];
s q38[1];
y q38[1];
cx q38[0],q38[1];
sdg q38[0];
h q38[0];
sdg q38[1];
h q38[1];
cx q38[0],q38[1];
h q38[0];
s q38[0];
h q38[1];
s q38[1];
cx q38[0],q38[1];
x q38[0];
h q38[0];
s q38[0];
cx q38[0],q38[1];
y q38[0];
h q38[0];
y q38[1];
h q38[1];
cx q38[0],q38[1];
h q38[0];
s q38[0];
sdg q38[1];
h q38[1];
cx q38[0],q38[1];
h q38[0];
s q38[0];
h q38[1];
s q38[1];
z q38[1];
s q38[1];
cx q38[0],q38[1];
h q38[0];
s q38[0];
h q38[1];
s q38[1];
cx q38[0],q38[1];
x q38[0];
s q38[0];
x q38[1];
cx q38[0],q38[1];
sdg q38[0];
h q38[0];
z q38[0];
sdg q38[0];
h q38[0];
h q38[1];
s q38[1];
z q38[1];
s q38[1];
h q38[1];
s q38[1];
cx q38[0],q38[1];
x q38[0];
s q38[0];
h q38[1];
s q38[1];
x q38[1];
s q38[1];
cx q38[0],q38[1];
x q38[0];
s q38[0];
sdg q38[1];
h q38[1];
z q38[1];
s q38[1];
h q38[1];
s q38[1];
cx q38[0],q38[1];
h q38[0];
s q38[0];
h q38[1];
s q38[1];
cx q38[0],q38[1];
sdg q38[0];
h q38[0];
z q38[1];
sdg q38[1];
h q38[1];
cx q38[0],q38[1];
sdg q38[0];
h q38[0];
sdg q38[1];
h q38[1];
cx q38[0],q38[1];
h q38[0];
s q38[0];
x q38[0];
h q38[0];
y q38[1];
s q38[1];
h q38[1];
s q38[1];
cx q38[0],q38[1];
y q38[0];
h q38[1];
s q38[1];
y q38[1];
h q38[1];
cx q38[0],q38[1];
h q38[0];
s q38[0];
sdg q38[1];
h q38[1];
cx q38[0],q38[1];
x q38[0];
s q38[0];
h q38[0];
s q38[0];
h q38[1];
s q38[1];
cx q38[0],q38[1];
h q38[0];
s q38[0];
h q38[1];
s q38[1];
cx q38[0],q38[1];
y q38[0];
s q38[0];
z q38[1];
h q38[1];
cx q38[0],q38[1];
h q38[0];
s q38[0];
y q38[0];
h q38[0];
s q38[1];
h q38[1];
s q38[1];
cx q38[0],q38[1];
sdg q38[0];
h q38[0];
h q38[1];
s q38[1];
cx q38[0],q38[1];
z q38[0];
h q38[0];
s q38[0];
z q38[1];
h q38[1];
s q38[1];
cx q38[0],q38[1];
h q38[0];
s q38[0];
h q38[1];
s q38[1];
cx q38[0],q38[1];
sdg q38[0];
h q38[0];
y q38[1];
h q38[1];
s q38[1];
cx q38[0],q38[1];
sdg q38[0];
h q38[0];
h q38[1];
s q38[1];
cx q38[0],q38[1];
h q38[0];
s q38[0];
sdg q38[0];
h q38[0];
h q38[1];
cx q38[0],q38[1];
sdg q38[0];
h q38[0];
y q38[0];
h q38[0];
y q38[1];
s q38[1];
cx q38[0],q38[1];
sdg q38[0];
h q38[0];
x q38[0];
s q38[0];
h q38[0];
s q38[0];
h q38[1];
s q38[1];
h q38[1];
cx q38[0],q38[1];
h q38[0];
s q38[0];
sdg q38[1];
h q38[1];
cx q38[0],q38[1];
z q38[0];
h q38[1];
s q38[1];
z q38[1];
h q38[1];
cx q38[0],q38[1];
sdg q38[0];
h q38[0];
sdg q38[1];
h q38[1];
cx q38[0],q38[1];
h q38[0];
s q38[0];
h q38[1];
s q38[1];
x q38[1];
s q38[1];
cx q38[0],q38[1];
z q38[0];
h q38[0];
s q38[0];
y q38[1];
h q38[1];
s q38[1];
cx q38[0],q38[1];
sdg q38[0];
h q38[0];
z q38[0];
sdg q38[1];
h q38[1];
x q38[1];
h q38[1];
s q38[1];
cx q38[0],q38[1];
sdg q38[0];
h q38[0];
sdg q38[1];
h q38[1];
cx q38[0],q38[1];
h q38[0];
s q38[0];
x q38[0];
sdg q38[0];
h q38[0];
x q38[1];
s q38[1];
h q38[1];
s q38[1];
cx q38[0],q38[1];
h q38[0];
s q38[0];
z q38[0];
s q38[0];
h q38[1];
s q38[1];
cx q38[0],q38[1];
sdg q38[0];
h q38[0];
z q38[1];
h q38[1];
s q38[1];
cx q38[0],q38[1];
h q38[0];
s q38[0];
h q38[1];
s q38[1];
cx q38[0],q38[1];
sdg q38[0];
h q38[0];
y q38[1];
cx q38[0],q38[1];
h q38[0];
s q38[0];
sdg q38[1];
h q38[1];
cx q38[0],q38[1];
x q38[0];
h q38[1];
s q38[1];
x q38[1];
h q38[1];
cx q38[0],q38[1];
sdg q38[0];
h q38[0];
h q38[1];
s q38[1];
cx q38[0],q38[1];
y q38[0];
s q38[0];
h q38[0];
s q38[0];
z q38[1];
h q38[1];
cx q38[0],q38[1];
sdg q38[0];
h q38[0];
sdg q38[1];
h q38[1];
cx q38[0],q38[1];
h q38[0];
s q38[0];
x q38[0];
s q38[0];
h q38[0];
s q38[0];
h q38[1];
s q38[1];
h q38[1];
s q38[1];
cx q38[0],q38[1];
sdg q38[0];
h q38[0];
z q38[0];
s q38[0];
y q38[1];
h q38[1];
s q38[1];
cx q38[0],q38[1];
z q38[0];
h q38[0];
s q38[0];
h q38[1];
s q38[1];
sdg q38[1];
h q38[1];
cx q38[0],q38[1];
sdg q38[0];
h q38[0];
x q38[0];
sdg q38[0];
h q38[0];
h q38[1];
s q38[1];
x q38[1];
cx q38[0],q38[1];
x q38[0];
sdg q38[0];
h q38[0];
sdg q38[1];
h q38[1];
y q38[1];
h q38[1];
s q38[1];
cx q38[0],q38[1];
sdg q38[0];
h q38[0];
sdg q38[1];
h q38[1];
cx q38[0],q38[1];
y q38[0];
h q38[0];
s q38[0];
x q38[1];
h q38[1];
cx q38[0],q38[1];
sdg q38[0];
h q38[0];
h q38[1];
s q38[1];
cx q38[0],q38[1];
h q38[0];
s q38[0];
s q38[0];
z q38[1];
h q38[1];
s q38[1];
cx q38[0],q38[1];
sdg q38[0];
h q38[0];
sdg q38[1];
h q38[1];
cx q38[0],q38[1];
y q38[0];
h q38[0];
y q38[1];
cx q38[0],q38[1];
sdg q38[0];
h q38[0];
z q38[0];
s q38[0];
h q38[0];
s q38[0];
sdg q38[1];
h q38[1];
s q38[1];
cx q38[0],q38[1];
sdg q38[0];
h q38[0];
sdg q38[1];
h q38[1];
cx q38[0],q38[1];
h q38[0];
s q38[0];
z q38[0];
h q38[0];
x q38[0];
h q38[0];
s q38[0];
z q38[1];
h q38[1];
s q38[1];
x q38[1];
h q38[1];
cx q38[0],q38[1];
sdg q38[0];
h q38[0];
sdg q38[1];
h q38[1];
cx q38[0],q38[1];
h q38[0];
s q38[0];
z q38[0];
s q38[0];
h q38[1];
s q38[1];
y q38[1];
sdg q38[1];
h q38[1];
cx q38[0],q38[1];
sdg q38[0];
h q38[0];
sdg q38[1];
h q38[1];
cx q38[0],q38[1];
h q38[0];
s q38[0];
z q38[0];
h q38[1];
s q38[1];
x q38[1];
h q38[1];
s q38[1];
cx q38[0],q38[1];
sdg q38[0];
h q38[0];
y q38[0];
h q38[0];
s q38[0];
sdg q38[1];
h q38[1];
x q38[1];
s q38[1];
cx q38[0],q38[1];
h q38[0];
s q38[0];
sdg q38[1];
h q38[1];
cx q38[0],q38[1];
x q38[0];
sdg q38[0];
h q38[0];
x q38[1];
h q38[1];
cx q38[0],q38[1];
h q38[0];
s q38[0];
z q38[0];
x q38[1];
s q38[1];
cx q38[0],q38[1];
sdg q38[0];
h q38[0];
h q38[1];
s q38[1];
cx q38[0],q38[1];
z q38[0];
z q38[1];
h q38[1];
s q38[1];
cx q38[0],q38[1];
h q38[0];
s q38[0];
y q38[0];
h q38[0];
s q38[0];
y q38[0];
x q38[1];
h q38[1];
s q38[1];
cx q38[0],q38[1];
y q38[0];
h q38[0];
s q38[0];
sdg q38[1];
h q38[1];
y q38[1];
h q38[1];
s q38[1];
cx q38[0],q38[1];
h q38[0];
s q38[0];
x q38[0];
sdg q38[0];
h q38[0];
h q38[1];
s q38[1];
h q38[1];
s q38[1];
cx q38[0],q38[1];
h q38[0];
s q38[0];
h q38[1];
s q38[1];
cx q38[0],q38[1];
z q38[0];
sdg q38[0];
h q38[0];
x q38[1];
sdg q38[1];
h q38[1];
cx q38[0],q38[1];
h q38[0];
s q38[0];
x q38[0];
s q38[0];
sdg q38[1];
h q38[1];
z q38[1];
s q38[1];
cx q38[0],q38[1];
sdg q38[0];
h q38[0];
z q38[0];
sdg q38[0];
h q38[0];
sdg q38[1];
h q38[1];
s q38[1];
h q38[1];
s q38[1];
cx q38[0],q38[1];
y q38[0];
h q38[0];
h q38[1];
s q38[1];
x q38[1];
cx q38[0],q38[1];
sdg q38[0];
h q38[0];
h q38[1];
s q38[1];
cx q38[0],q38[1];
h q38[0];
s q38[0];
z q38[0];
s q38[0];
y q38[1];
sdg q38[1];
h q38[1];
cx q38[0],q38[1];
h q38[0];
s q38[0];
sdg q38[1];
h q38[1];
cx q38[0],q38[1];
x q38[0];
s q38[0];
h q38[1];
s q38[1];
z q38[1];
h q38[1];
s q38[1];
cx q38[0],q38[1];
x q38[0];
s q38[0];
h q38[0];
s q38[0];
h q38[1];
s q38[1];
s q38[1];
h q38[1];
s q38[1];
cx q38[0],q38[1];
sdg q38[0];
h q38[0];
sdg q38[1];
h q38[1];
cx q38[0],q38[1];
h q38[0];
s q38[0];
y q38[0];
sdg q38[0];
h q38[0];
h q38[1];
s q38[1];
cx q38[0],q38[1];
x q38[0];
s q38[0];
y q38[1];
h q38[1];
cx q38[0],q38[1];
h q38[0];
h q38[1];
s q38[1];
x q38[1];
h q38[1];
cx q38[0],q38[1];
h q38[0];
s q38[0];
z q38[0];
s q38[0];
h q38[0];
s q38[0];
sdg q38[1];
h q38[1];
x q38[1];
s q38[1];
cx q38[0],q38[1];
sdg q38[0];
h q38[0];
sdg q38[1];
h q38[1];
cx q38[0],q38[1];
h q38[0];
s q38[0];
h q38[1];
s q38[1];
cx q38[0],q38[1];
h q38[0];
y q38[1];
s q38[1];
cx q38[0],q38[1];
sdg q38[0];
h q38[0];
h q38[1];
s q38[1];
y q38[1];
cx q38[0],q38[1];
sdg q38[0];
h q38[0];
h q38[1];
s q38[1];
cx q38[0],q38[1];
h q38[0];
s q38[0];
y q38[0];
h q38[0];
s q38[0];
cx q38[0],q38[1];
h q38[0];
s q38[0];
h q38[1];
s q38[1];
cx q38[0],q38[1];
y q38[0];
h q38[0];
x q38[1];
s q38[1];
h q38[1];
s q38[1];
cx q38[0],q38[1];
h q38[0];
s q38[0];
sdg q38[1];
h q38[1];
cx q38[0],q38[1];
y q38[0];
h q38[0];
s q38[0];
h q38[1];
s q38[1];
z q38[1];
s q38[1];
h q38[1];
s q38[1];
cx q38[0],q38[1];
sdg q38[0];
h q38[0];
sdg q38[1];
h q38[1];
cx q38[0],q38[1];
h q38[0];
s q38[0];
s q38[0];
h q38[0];
s q38[0];
y q38[1];
s q38[1];
h q38[1];
s q38[1];
cx q38[0],q38[1];
h q38[0];
s q38[0];
y q38[0];
z q38[1];
s q38[1];
cx q38[0],q38[1];
h q38[0];
s q38[0];
y q38[0];
s q38[0];
h q38[0];
s q38[0];
h q38[1];
s q38[1];
sdg q38[1];
h q38[1];
cx q38[0],q38[1];
sdg q38[0];
h q38[0];
sdg q38[1];
h q38[1];
cx q38[0],q38[1];
z q38[0];
h q38[0];
h q38[1];
s q38[1];
z q38[1];
s q38[1];
h q38[1];
s q38[1];
cx q38[0],q38[1];
sdg q38[0];
h q38[0];
s q38[0];
z q38[1];
cx q38[0],q38[1];
sdg q38[0];
h q38[0];
z q38[0];
y q38[0];
h q38[0];
y q38[1];
h q38[1];
s q38[1];
y q38[1];
s q38[1];
h q38[1];
s q38[1];
cx q38[0],q38[1];
h q38[0];
s q38[0];
h q38[1];
s q38[1];
cx q38[0],q38[1];
h q38[0];
s q38[0];
y q38[0];
s q38[0];
h q38[0];
s q38[0];
z q38[1];
s q38[1];
cx q38[0],q38[1];
z q38[0];
sdg q38[0];
h q38[0];
z q38[1];
s q38[1];
cx q38[0],q38[1];
y q38[0];
s q38[0];
h q38[0];
s q38[0];
sdg q38[1];
h q38[1];
x q38[1];
cx q38[0],q38[1];
h q38[0];
s q38[0];
sdg q38[1];
h q38[1];
cx q38[0],q38[1];
sdg q38[0];
h q38[0];
h q38[1];
s q38[1];
y q38[1];
s q38[1];
h q38[1];
s q38[1];
cx q38[0],q38[1];
s q38[0];
h q38[0];
s q38[0];
h q38[1];
s q38[1];
x q38[1];
sdg q38[1];
h q38[1];
cx q38[0],q38[1];
sdg q38[0];
h q38[0];
h q38[1];
s q38[1];
cx q38[0],q38[1];
h q38[0];
s q38[0];
x q38[0];
s q38[0];
y q38[1];
cx q38[0],q38[1];
sdg q38[0];
h q38[0];
sdg q38[1];
h q38[1];
cx q38[0],q38[1];
h q38[0];
s q38[0];
h q38[1];
s q38[1];
cx q38[0],q38[1];
y q38[0];
h q38[0];
s q38[0];
x q38[1];
h q38[1];
s q38[1];
cx q38[0],q38[1];
h q38[0];
s q38[0];
sdg q38[1];
h q38[1];
cx q38[0],q38[1];
