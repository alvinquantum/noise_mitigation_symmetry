OPENQASM 2.0;
include "qelib1.inc";
qreg q26[2];
y q26[0];
s q26[0];
h q26[0];
s q26[0];
sdg q26[1];
h q26[1];
cx q26[0],q26[1];
sdg q26[0];
h q26[0];
h q26[1];
s q26[1];
cx q26[0],q26[1];
y q26[0];
h q26[0];
s q26[1];
cx q26[0],q26[1];
sdg q26[0];
h q26[0];
s q26[0];
h q26[0];
s q26[0];
y q26[1];
sdg q26[1];
h q26[1];
cx q26[0],q26[1];
sdg q26[0];
h q26[0];
sdg q26[1];
h q26[1];
cx q26[0],q26[1];
x q26[0];
x q26[1];
s q26[1];
cx q26[0],q26[1];
sdg q26[0];
h q26[0];
sdg q26[1];
h q26[1];
cx q26[0],q26[1];
h q26[0];
s q26[0];
y q26[0];
sdg q26[0];
h q26[0];
y q26[1];
cx q26[0],q26[1];
x q26[0];
s q26[0];
sdg q26[1];
h q26[1];
y q26[1];
s q26[1];
cx q26[0],q26[1];
h q26[0];
s q26[0];
h q26[0];
sdg q26[1];
h q26[1];
z q26[1];
sdg q26[1];
h q26[1];
cx q26[0],q26[1];
sdg q26[0];
h q26[0];
sdg q26[1];
h q26[1];
cx q26[0],q26[1];
y q26[0];
s q26[0];
z q26[1];
cx q26[0],q26[1];
sdg q26[0];
h q26[0];
z q26[0];
s q26[0];
sdg q26[1];
h q26[1];
z q26[1];
s q26[1];
h q26[1];
s q26[1];
cx q26[0],q26[1];
sdg q26[0];
h q26[0];
sdg q26[1];
h q26[1];
cx q26[0],q26[1];
h q26[0];
s q26[0];
z q26[0];
h q26[0];
s q26[0];
x q26[1];
sdg q26[1];
h q26[1];
cx q26[0],q26[1];
h q26[0];
s q26[0];
sdg q26[1];
h q26[1];
cx q26[0],q26[1];
x q26[0];
h q26[0];
y q26[1];
sdg q26[1];
h q26[1];
cx q26[0],q26[1];
sdg q26[0];
h q26[0];
y q26[0];
sdg q26[0];
h q26[0];
h q26[1];
s q26[1];
z q26[1];
s q26[1];
h q26[1];
s q26[1];
cx q26[0],q26[1];
sdg q26[0];
h q26[0];
sdg q26[1];
h q26[1];
cx q26[0],q26[1];
h q26[0];
s q26[0];
h q26[1];
s q26[1];
cx q26[0],q26[1];
x q26[0];
s q26[0];
cx q26[0],q26[1];
h q26[0];
s q26[0];
z q26[0];
h q26[0];
s q26[0];
h q26[1];
s q26[1];
z q26[1];
s q26[1];
h q26[1];
s q26[1];
cx q26[0],q26[1];
y q26[0];
s q26[0];
h q26[0];
s q26[0];
sdg q26[1];
h q26[1];
y q26[1];
cx q26[0],q26[1];
sdg q26[0];
h q26[0];
h q26[1];
s q26[1];
cx q26[0],q26[1];
h q26[0];
s q26[0];
x q26[1];
s q26[1];
h q26[1];
s q26[1];
cx q26[0],q26[1];
sdg q26[0];
h q26[0];
h q26[1];
s q26[1];
cx q26[0],q26[1];
z q26[0];
sdg q26[0];
h q26[0];
cx q26[0],q26[1];
h q26[0];
s q26[0];
sdg q26[1];
h q26[1];
cx q26[0],q26[1];
h q26[0];
s q26[0];
s q26[1];
h q26[1];
s q26[1];
cx q26[0],q26[1];
h q26[0];
s q26[0];
h q26[1];
s q26[1];
cx q26[0],q26[1];
y q26[0];
s q26[0];
sdg q26[0];
h q26[0];
h q26[1];
y q26[1];
h q26[1];
cx q26[0],q26[1];
sdg q26[0];
h q26[0];
sdg q26[1];
h q26[1];
cx q26[0],q26[1];
x q26[0];
h q26[0];
s q26[0];
h q26[1];
s q26[1];
cx q26[0],q26[1];
sdg q26[0];
h q26[0];
z q26[0];
s q26[0];
h q26[1];
s q26[1];
sdg q26[1];
h q26[1];
cx q26[0],q26[1];
sdg q26[0];
h q26[0];
sdg q26[1];
h q26[1];
cx q26[0],q26[1];
h q26[0];
s q26[0];
z q26[0];
h q26[0];
s q26[0];
x q26[1];
h q26[1];
s q26[1];
cx q26[0],q26[1];
sdg q26[0];
h q26[0];
sdg q26[1];
h q26[1];
cx q26[0],q26[1];
s q26[0];
h q26[0];
s q26[0];
h q26[1];
s q26[1];
h q26[1];
s q26[1];
cx q26[0],q26[1];
h q26[0];
s q26[0];
y q26[0];
h q26[0];
s q26[0];
sdg q26[1];
h q26[1];
sdg q26[1];
h q26[1];
cx q26[0],q26[1];
sdg q26[0];
h q26[0];
sdg q26[1];
h q26[1];
cx q26[0],q26[1];
h q26[0];
s q26[0];
y q26[0];
h q26[0];
s q26[0];
h q26[1];
s q26[1];
x q26[1];
h q26[1];
cx q26[0],q26[1];
h q26[0];
s q26[0];
h q26[1];
s q26[1];
cx q26[0],q26[1];
y q26[0];
sdg q26[0];
h q26[0];
s q26[1];
h q26[1];
s q26[1];
cx q26[0],q26[1];
sdg q26[0];
h q26[0];
sdg q26[1];
h q26[1];
cx q26[0],q26[1];
h q26[0];
s q26[0];
h q26[1];
s q26[1];
z q26[1];
sdg q26[1];
h q26[1];
cx q26[0],q26[1];
sdg q26[0];
h q26[0];
sdg q26[1];
h q26[1];
cx q26[0],q26[1];
y q26[0];
s q26[0];
h q26[1];
s q26[1];
y q26[1];
s q26[1];
cx q26[0],q26[1];
h q26[0];
s q26[0];
sdg q26[1];
h q26[1];
z q26[1];
s q26[1];
h q26[1];
s q26[1];
cx q26[0],q26[1];
h q26[0];
s q26[0];
h q26[1];
s q26[1];
cx q26[0],q26[1];
z q26[0];
s q26[1];
h q26[1];
s q26[1];
cx q26[0],q26[1];
sdg q26[0];
h q26[0];
sdg q26[1];
h q26[1];
cx q26[0],q26[1];
z q26[0];
sdg q26[0];
h q26[0];
z q26[1];
h q26[1];
s q26[1];
cx q26[0],q26[1];
h q26[0];
s q26[0];
z q26[0];
s q26[0];
h q26[1];
s q26[1];
x q26[1];
s q26[1];
cx q26[0],q26[1];
sdg q26[0];
h q26[0];
h q26[1];
s q26[1];
cx q26[0],q26[1];
h q26[0];
s q26[0];
y q26[0];
z q26[1];
h q26[1];
s q26[1];
cx q26[0],q26[1];
h q26[0];
s q26[0];
sdg q26[1];
h q26[1];
cx q26[0],q26[1];
y q26[0];
sdg q26[0];
h q26[0];
h q26[1];
s q26[1];
x q26[1];
h q26[1];
s q26[1];
cx q26[0],q26[1];
sdg q26[0];
h q26[0];
z q26[0];
s q26[0];
h q26[0];
s q26[0];
sdg q26[1];
h q26[1];
z q26[1];
sdg q26[1];
h q26[1];
cx q26[0],q26[1];
h q26[0];
s q26[0];
h q26[1];
s q26[1];
cx q26[0],q26[1];
x q26[1];
h q26[1];
cx q26[0],q26[1];
h q26[0];
s q26[0];
x q26[0];
s q26[0];
h q26[0];
s q26[0];
y q26[1];
sdg q26[1];
h q26[1];
cx q26[0],q26[1];
z q26[0];
sdg q26[1];
h q26[1];
s q26[1];
cx q26[0],q26[1];
z q26[0];
h q26[0];
z q26[1];
h q26[1];
s q26[1];
cx q26[0],q26[1];
h q26[0];
s q26[0];
sdg q26[1];
h q26[1];
cx q26[0],q26[1];
x q26[0];
s q26[0];
x q26[1];
sdg q26[1];
h q26[1];
cx q26[0],q26[1];
h q26[0];
s q26[0];
h q26[1];
s q26[1];
cx q26[0],q26[1];
sdg q26[0];
h q26[0];
z q26[1];
s q26[1];
cx q26[0],q26[1];
sdg q26[0];
h q26[0];
x q26[0];
sdg q26[0];
h q26[0];
x q26[1];
s q26[1];
cx q26[0],q26[1];
sdg q26[0];
h q26[0];
z q26[0];
h q26[1];
s q26[1];
z q26[1];
sdg q26[1];
h q26[1];
cx q26[0],q26[1];
h q26[0];
s q26[0];
x q26[0];
h q26[1];
s q26[1];
x q26[1];
h q26[1];
cx q26[0],q26[1];
h q26[0];
s q26[0];
y q26[0];
sdg q26[0];
h q26[0];
z q26[1];
cx q26[0],q26[1];
sdg q26[0];
h q26[0];
sdg q26[1];
h q26[1];
cx q26[0],q26[1];
s q26[0];
y q26[1];
s q26[1];
h q26[1];
s q26[1];
cx q26[0],q26[1];
h q26[0];
s q26[0];
sdg q26[1];
h q26[1];
cx q26[0],q26[1];
x q26[0];
sdg q26[0];
h q26[0];
y q26[1];
s q26[1];
h q26[1];
s q26[1];
cx q26[0],q26[1];
h q26[0];
s q26[0];
z q26[0];
s q26[0];
sdg q26[1];
h q26[1];
s q26[1];
h q26[1];
s q26[1];
cx q26[0],q26[1];
sdg q26[0];
h q26[0];
sdg q26[1];
h q26[1];
cx q26[0],q26[1];
h q26[0];
s q26[0];
x q26[0];
s q26[0];
h q26[0];
s q26[0];
y q26[1];
s q26[1];
h q26[1];
s q26[1];
cx q26[0],q26[1];
y q26[0];
s q26[0];
h q26[1];
s q26[1];
x q26[1];
h q26[1];
s q26[1];
cx q26[0],q26[1];
sdg q26[0];
h q26[0];
sdg q26[1];
h q26[1];
cx q26[0],q26[1];
h q26[0];
s q26[0];
y q26[0];
s q26[0];
h q26[0];
s q26[0];
h q26[1];
s q26[1];
y q26[1];
s q26[1];
cx q26[0],q26[1];
sdg q26[0];
h q26[0];
z q26[0];
h q26[0];
s q26[0];
s q26[1];
h q26[1];
s q26[1];
cx q26[0],q26[1];
sdg q26[0];
h q26[0];
sdg q26[1];
h q26[1];
cx q26[0],q26[1];
y q26[0];
s q26[0];
h q26[0];
s q26[0];
h q26[1];
s q26[1];
y q26[1];
s q26[1];
h q26[1];
s q26[1];
cx q26[0],q26[1];
x q26[0];
sdg q26[0];
h q26[0];
sdg q26[1];
h q26[1];
y q26[1];
h q26[1];
s q26[1];
cx q26[0],q26[1];
sdg q26[0];
h q26[0];
y q26[1];
cx q26[0],q26[1];
sdg q26[0];
h q26[0];
sdg q26[1];
h q26[1];
cx q26[0],q26[1];
sdg q26[0];
h q26[0];
h q26[1];
s q26[1];
y q26[1];
h q26[1];
cx q26[0],q26[1];
sdg q26[0];
h q26[0];
sdg q26[1];
h q26[1];
cx q26[0],q26[1];
y q26[0];
s q26[0];
h q26[0];
s q26[0];
h q26[1];
s q26[1];
sdg q26[1];
h q26[1];
cx q26[0],q26[1];
h q26[0];
s q26[0];
x q26[0];
s q26[0];
h q26[0];
s q26[0];
h q26[1];
s q26[1];
y q26[1];
s q26[1];
cx q26[0],q26[1];
x q26[0];
s q26[0];
s q26[0];
h q26[1];
s q26[1];
y q26[1];
h q26[1];
s q26[1];
h q26[1];
s q26[1];
cx q26[0],q26[1];
sdg q26[0];
h q26[0];
sdg q26[1];
h q26[1];
cx q26[0],q26[1];
s q26[0];
x q26[1];
h q26[1];
cx q26[0],q26[1];
h q26[0];
s q26[0];
y q26[0];
h q26[0];
s q26[0];
sdg q26[1];
h q26[1];
x q26[1];
cx q26[0],q26[1];
sdg q26[0];
h q26[0];
h q26[1];
s q26[1];
cx q26[0],q26[1];
h q26[0];
s q26[0];
y q26[0];
sdg q26[0];
h q26[0];
z q26[1];
sdg q26[1];
h q26[1];
cx q26[0],q26[1];
h q26[0];
s q26[0];
y q26[0];
s q26[0];
h q26[0];
s q26[0];
h q26[1];
s q26[1];
y q26[1];
h q26[1];
s q26[1];
cx q26[0],q26[1];
h q26[0];
s q26[0];
sdg q26[1];
h q26[1];
cx q26[0],q26[1];
h q26[0];
x q26[1];
sdg q26[1];
h q26[1];
cx q26[0],q26[1];
sdg q26[0];
h q26[0];
z q26[0];
h q26[0];
s q26[0];
h q26[1];
s q26[1];
sdg q26[1];
h q26[1];
cx q26[0],q26[1];
h q26[0];
s q26[0];
h q26[1];
s q26[1];
cx q26[0],q26[1];
s q26[0];
h q26[0];
s q26[0];
z q26[0];
sdg q26[0];
h q26[0];
x q26[1];
h q26[1];
s q26[1];
y q26[1];
s q26[1];
h q26[1];
s q26[1];
cx q26[0],q26[1];
sdg q26[0];
h q26[0];
h q26[1];
s q26[1];
cx q26[0],q26[1];
h q26[0];
x q26[0];
h q26[0];
z q26[1];
sdg q26[1];
h q26[1];
x q26[1];
h q26[1];
cx q26[0],q26[1];
h q26[0];
s q26[0];
y q26[0];
s q26[0];
h q26[0];
s q26[0];
sdg q26[1];
h q26[1];
z q26[1];
h q26[1];
s q26[1];
cx q26[0],q26[1];
sdg q26[0];
h q26[0];
h q26[1];
s q26[1];
cx q26[0],q26[1];
z q26[0];
s q26[1];
cx q26[0],q26[1];
sdg q26[0];
h q26[0];
x q26[0];
x q26[1];
s q26[1];
cx q26[0],q26[1];
s q26[0];
h q26[0];
s q26[0];
sdg q26[1];
h q26[1];
cx q26[0],q26[1];
sdg q26[0];
h q26[0];
sdg q26[1];
h q26[1];
cx q26[0],q26[1];
y q26[0];
s q26[0];
h q26[0];
s q26[0];
h q26[1];
s q26[1];
z q26[1];
sdg q26[1];
h q26[1];
cx q26[0],q26[1];
y q26[0];
h q26[0];
s q26[0];
z q26[1];
h q26[1];
cx q26[0],q26[1];
sdg q26[0];
h q26[0];
x q26[0];
sdg q26[0];
h q26[0];
sdg q26[1];
h q26[1];
z q26[1];
s q26[1];
h q26[1];
s q26[1];
cx q26[0],q26[1];
sdg q26[0];
h q26[0];
sdg q26[1];
h q26[1];
cx q26[0],q26[1];
h q26[0];
s q26[0];
h q26[1];
s q26[1];
y q26[1];
cx q26[0],q26[1];
x q26[0];
h q26[0];
h q26[1];
s q26[1];
x q26[1];
h q26[1];
cx q26[0],q26[1];
h q26[0];
s q26[0];
z q26[0];
h q26[0];
s q26[0];
z q26[1];
cx q26[0],q26[1];
h q26[0];
s q26[0];
y q26[0];
s q26[0];
h q26[0];
s q26[0];
y q26[1];
h q26[1];
cx q26[0],q26[1];
h q26[0];
s q26[0];
y q26[0];
h q26[1];
s q26[1];
cx q26[0],q26[1];
h q26[0];
s q26[0];
sdg q26[1];
h q26[1];
cx q26[0],q26[1];
s q26[0];
z q26[1];
sdg q26[1];
h q26[1];
cx q26[0],q26[1];
h q26[0];
s q26[0];
h q26[1];
s q26[1];
cx q26[0],q26[1];
x q26[0];
h q26[0];
s q26[0];
y q26[1];
s q26[1];
h q26[1];
s q26[1];
cx q26[0],q26[1];
z q26[0];
z q26[1];
h q26[1];
x q26[1];
h q26[1];
cx q26[0],q26[1];
sdg q26[0];
h q26[0];
sdg q26[1];
h q26[1];
cx q26[0],q26[1];
y q26[0];
h q26[0];
s q26[0];
z q26[1];
h q26[1];
cx q26[0],q26[1];
sdg q26[0];
h q26[0];
sdg q26[1];
h q26[1];
cx q26[0],q26[1];
h q26[1];
s q26[1];
