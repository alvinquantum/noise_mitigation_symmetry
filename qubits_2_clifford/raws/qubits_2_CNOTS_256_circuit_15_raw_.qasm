OPENQASM 2.0;
include "qelib1.inc";
qreg q14[2];
x q14[0];
z q14[0];
h q14[0];
h q14[1];
s q14[1];
s q14[1];
h q14[1];
s q14[1];
z q14[1];
s q14[1];
cx q14[0],q14[1];
sdg q14[0];
h q14[0];
h q14[1];
s q14[1];
cx q14[0],q14[1];
y q14[0];
s q14[0];
sdg q14[1];
h q14[1];
cx q14[0],q14[1];
sdg q14[0];
h q14[0];
z q14[0];
s q14[0];
z q14[1];
h q14[1];
cx q14[0],q14[1];
sdg q14[0];
h q14[0];
y q14[0];
h q14[0];
s q14[0];
y q14[1];
s q14[1];
h q14[1];
s q14[1];
cx q14[0],q14[1];
h q14[0];
s q14[0];
s q14[0];
h q14[0];
s q14[0];
sdg q14[1];
h q14[1];
z q14[1];
s q14[1];
cx q14[0],q14[1];
sdg q14[0];
h q14[0];
h q14[1];
s q14[1];
cx q14[0],q14[1];
h q14[0];
s q14[0];
y q14[0];
s q14[0];
h q14[0];
s q14[0];
s q14[1];
h q14[1];
s q14[1];
cx q14[0],q14[1];
sdg q14[0];
h q14[0];
h q14[1];
s q14[1];
cx q14[0],q14[1];
h q14[0];
s q14[0];
y q14[0];
s q14[0];
h q14[0];
s q14[0];
x q14[1];
cx q14[0],q14[1];
h q14[0];
s q14[0];
sdg q14[1];
h q14[1];
cx q14[0],q14[1];
z q14[0];
s q14[0];
cx q14[0],q14[1];
sdg q14[0];
h q14[0];
h q14[1];
s q14[1];
cx q14[0],q14[1];
h q14[0];
s q14[0];
y q14[0];
h q14[0];
s q14[0];
x q14[1];
h q14[1];
s q14[1];
cx q14[0],q14[1];
h q14[0];
s q14[0];
sdg q14[1];
h q14[1];
cx q14[0],q14[1];
s q14[0];
y q14[1];
s q14[1];
cx q14[0],q14[1];
sdg q14[0];
h q14[0];
z q14[0];
sdg q14[0];
h q14[0];
h q14[1];
s q14[1];
y q14[1];
h q14[1];
cx q14[0],q14[1];
sdg q14[0];
h q14[0];
sdg q14[1];
h q14[1];
cx q14[0],q14[1];
h q14[0];
s q14[0];
y q14[0];
s q14[0];
cx q14[0],q14[1];
sdg q14[0];
h q14[0];
sdg q14[1];
h q14[1];
cx q14[0],q14[1];
x q14[0];
sdg q14[0];
h q14[0];
h q14[0];
s q14[0];
h q14[1];
s q14[1];
z q14[1];
h q14[1];
s q14[1];
y q14[1];
s q14[1];
h q14[1];
s q14[1];
cx q14[0],q14[1];
sdg q14[0];
h q14[0];
h q14[1];
s q14[1];
cx q14[0],q14[1];
h q14[0];
s q14[0];
h q14[0];
z q14[1];
cx q14[0],q14[1];
h q14[0];
s q14[0];
x q14[1];
cx q14[0],q14[1];
h q14[0];
s q14[0];
sdg q14[1];
h q14[1];
cx q14[0],q14[1];
x q14[0];
h q14[1];
s q14[1];
y q14[1];
cx q14[0],q14[1];
sdg q14[0];
h q14[0];
sdg q14[1];
h q14[1];
cx q14[0],q14[1];
x q14[0];
s q14[0];
y q14[1];
sdg q14[1];
h q14[1];
cx q14[0],q14[1];
sdg q14[0];
h q14[0];
sdg q14[1];
h q14[1];
cx q14[0],q14[1];
h q14[0];
s q14[0];
y q14[0];
y q14[0];
h q14[1];
s q14[1];
y q14[1];
h q14[1];
s q14[1];
x q14[1];
h q14[1];
cx q14[0],q14[1];
sdg q14[0];
h q14[0];
sdg q14[1];
h q14[1];
cx q14[0],q14[1];
h q14[0];
s q14[0];
s q14[0];
h q14[1];
s q14[1];
y q14[1];
h q14[1];
cx q14[0],q14[1];
h q14[0];
s q14[0];
y q14[0];
s q14[0];
z q14[1];
sdg q14[1];
h q14[1];
cx q14[0],q14[1];
sdg q14[0];
h q14[0];
sdg q14[1];
h q14[1];
cx q14[0],q14[1];
h q14[0];
s q14[0];
h q14[1];
s q14[1];
cx q14[0],q14[1];
x q14[0];
s q14[0];
h q14[0];
s q14[0];
cx q14[0],q14[1];
sdg q14[0];
h q14[0];
sdg q14[1];
h q14[1];
cx q14[0],q14[1];
h q14[0];
s q14[0];
z q14[0];
h q14[1];
s q14[1];
x q14[1];
sdg q14[1];
h q14[1];
cx q14[0],q14[1];
sdg q14[1];
h q14[1];
s q14[1];
cx q14[0],q14[1];
sdg q14[0];
h q14[0];
sdg q14[1];
h q14[1];
cx q14[0],q14[1];
x q14[0];
sdg q14[1];
h q14[1];
cx q14[0],q14[1];
sdg q14[0];
h q14[0];
h q14[1];
s q14[1];
cx q14[0],q14[1];
h q14[0];
s q14[0];
y q14[0];
s q14[0];
y q14[0];
sdg q14[0];
h q14[0];
sdg q14[1];
h q14[1];
x q14[1];
s q14[1];
h q14[1];
s q14[1];
cx q14[0],q14[1];
h q14[0];
s q14[0];
z q14[0];
s q14[0];
h q14[0];
s q14[0];
x q14[0];
h q14[0];
s q14[0];
h q14[1];
s q14[1];
z q14[1];
s q14[1];
sdg q14[1];
h q14[1];
cx q14[0],q14[1];
x q14[0];
h q14[0];
s q14[0];
y q14[1];
sdg q14[1];
h q14[1];
cx q14[0],q14[1];
sdg q14[0];
h q14[0];
sdg q14[1];
h q14[1];
cx q14[0],q14[1];
h q14[0];
s q14[0];
z q14[0];
s q14[0];
h q14[0];
s q14[0];
x q14[1];
s q14[1];
h q14[1];
s q14[1];
cx q14[0],q14[1];
sdg q14[0];
h q14[0];
sdg q14[1];
h q14[1];
cx q14[0],q14[1];
h q14[0];
s q14[0];
h q14[1];
s q14[1];
cx q14[0],q14[1];
h q14[0];
s q14[0];
x q14[1];
h q14[1];
s q14[1];
cx q14[0],q14[1];
h q14[0];
s q14[0];
h q14[1];
s q14[1];
z q14[1];
cx q14[0],q14[1];
z q14[0];
h q14[1];
s q14[1];
x q14[1];
h q14[1];
cx q14[0],q14[1];
h q14[0];
s q14[0];
sdg q14[1];
h q14[1];
cx q14[0],q14[1];
z q14[0];
h q14[0];
h q14[1];
s q14[1];
z q14[1];
cx q14[0],q14[1];
sdg q14[0];
h q14[0];
h q14[1];
s q14[1];
cx q14[0],q14[1];
h q14[0];
s q14[0];
y q14[0];
s q14[0];
h q14[0];
s q14[0];
x q14[1];
s q14[1];
cx q14[0],q14[1];
h q14[0];
s q14[0];
sdg q14[1];
h q14[1];
cx q14[0],q14[1];
sdg q14[0];
h q14[0];
y q14[1];
cx q14[0],q14[1];
sdg q14[0];
h q14[0];
x q14[0];
sdg q14[1];
h q14[1];
x q14[1];
cx q14[0],q14[1];
h q14[0];
s q14[0];
y q14[0];
sdg q14[1];
h q14[1];
x q14[1];
h q14[1];
s q14[1];
cx q14[0],q14[1];
h q14[0];
s q14[0];
y q14[0];
h q14[1];
s q14[1];
z q14[1];
sdg q14[1];
h q14[1];
cx q14[0],q14[1];
sdg q14[0];
h q14[0];
y q14[0];
sdg q14[0];
h q14[0];
sdg q14[1];
h q14[1];
z q14[1];
s q14[1];
h q14[1];
s q14[1];
cx q14[0],q14[1];
sdg q14[0];
h q14[0];
z q14[0];
h q14[0];
s q14[0];
s q14[1];
cx q14[0],q14[1];
sdg q14[0];
h q14[0];
sdg q14[1];
h q14[1];
cx q14[0],q14[1];
x q14[0];
sdg q14[0];
h q14[0];
h q14[1];
s q14[1];
z q14[1];
cx q14[0],q14[1];
z q14[0];
h q14[0];
s q14[0];
x q14[1];
s q14[1];
h q14[1];
s q14[1];
cx q14[0],q14[1];
h q14[0];
s q14[0];
sdg q14[0];
h q14[0];
h q14[1];
s q14[1];
x q14[1];
cx q14[0],q14[1];
y q14[0];
h q14[0];
s q14[0];
h q14[1];
s q14[1];
z q14[1];
h q14[1];
cx q14[0],q14[1];
sdg q14[0];
h q14[0];
sdg q14[1];
h q14[1];
cx q14[0],q14[1];
x q14[0];
s q14[0];
h q14[0];
s q14[0];
x q14[1];
cx q14[0],q14[1];
z q14[0];
h q14[0];
h q14[1];
s q14[1];
cx q14[0],q14[1];
z q14[0];
s q14[0];
h q14[0];
s q14[0];
h q14[1];
cx q14[0],q14[1];
h q14[0];
s q14[0];
h q14[1];
s q14[1];
cx q14[0],q14[1];
x q14[0];
h q14[0];
s q14[0];
y q14[1];
s q14[1];
h q14[1];
s q14[1];
cx q14[0],q14[1];
z q14[0];
s q14[0];
h q14[0];
s q14[0];
x q14[1];
s q14[1];
h q14[1];
s q14[1];
cx q14[0],q14[1];
sdg q14[0];
h q14[0];
sdg q14[1];
h q14[1];
cx q14[0],q14[1];
z q14[0];
h q14[0];
y q14[1];
s q14[1];
cx q14[0],q14[1];
h q14[0];
s q14[0];
h q14[0];
s q14[0];
h q14[1];
s q14[1];
h q14[1];
s q14[1];
cx q14[0],q14[1];
sdg q14[0];
h q14[0];
y q14[0];
h q14[0];
s q14[0];
z q14[1];
h q14[1];
s q14[1];
cx q14[0],q14[1];
sdg q14[0];
h q14[0];
sdg q14[1];
h q14[1];
cx q14[0],q14[1];
z q14[0];
sdg q14[0];
h q14[0];
h q14[1];
s q14[1];
x q14[1];
sdg q14[1];
h q14[1];
cx q14[0],q14[1];
y q14[0];
s q14[0];
h q14[0];
s q14[0];
sdg q14[1];
h q14[1];
cx q14[0],q14[1];
z q14[0];
h q14[0];
h q14[1];
s q14[1];
y q14[1];
sdg q14[1];
h q14[1];
cx q14[0],q14[1];
sdg q14[0];
h q14[0];
sdg q14[1];
h q14[1];
cx q14[0],q14[1];
h q14[0];
s q14[0];
x q14[0];
h q14[0];
h q14[1];
s q14[1];
s q14[1];
cx q14[0],q14[1];
h q14[0];
s q14[0];
h q14[1];
s q14[1];
cx q14[0],q14[1];
x q14[0];
s q14[0];
x q14[1];
s q14[1];
h q14[1];
s q14[1];
cx q14[0],q14[1];
sdg q14[0];
h q14[0];
h q14[1];
s q14[1];
cx q14[0],q14[1];
x q14[0];
s q14[0];
h q14[0];
s q14[0];
x q14[1];
s q14[1];
cx q14[0],q14[1];
h q14[0];
s q14[0];
x q14[0];
h q14[0];
s q14[0];
y q14[0];
h q14[0];
x q14[1];
x q14[1];
sdg q14[1];
h q14[1];
cx q14[0],q14[1];
z q14[0];
s q14[0];
sdg q14[1];
h q14[1];
cx q14[0],q14[1];
z q14[0];
sdg q14[1];
h q14[1];
y q14[1];
cx q14[0],q14[1];
sdg q14[0];
h q14[0];
sdg q14[1];
h q14[1];
cx q14[0],q14[1];
h q14[0];
s q14[0];
z q14[0];
h q14[1];
cx q14[0],q14[1];
h q14[0];
s q14[0];
h q14[1];
s q14[1];
h q14[1];
s q14[1];
cx q14[0],q14[1];
sdg q14[0];
h q14[0];
sdg q14[1];
h q14[1];
cx q14[0],q14[1];
h q14[0];
s q14[0];
h q14[1];
s q14[1];
cx q14[0],q14[1];
z q14[0];
sdg q14[0];
h q14[0];
x q14[1];
sdg q14[1];
h q14[1];
cx q14[0],q14[1];
h q14[0];
s q14[0];
x q14[0];
s q14[0];
z q14[1];
h q14[1];
s q14[1];
cx q14[0],q14[1];
sdg q14[0];
h q14[0];
sdg q14[1];
h q14[1];
cx q14[0],q14[1];
h q14[1];
s q14[1];
x q14[1];
sdg q14[1];
h q14[1];
cx q14[0],q14[1];
h q14[0];
s q14[0];
z q14[0];
s q14[0];
h q14[0];
s q14[0];
h q14[1];
s q14[1];
y q14[1];
s q14[1];
cx q14[0],q14[1];
sdg q14[0];
h q14[0];
h q14[1];
s q14[1];
cx q14[0],q14[1];
y q14[0];
h q14[0];
x q14[1];
h q14[1];
s q14[1];
cx q14[0],q14[1];
s q14[0];
y q14[1];
s q14[1];
cx q14[0],q14[1];
z q14[1];
s q14[1];
cx q14[0],q14[1];
h q14[0];
s q14[0];
sdg q14[1];
h q14[1];
cx q14[0],q14[1];
y q14[0];
h q14[0];
sdg q14[1];
h q14[1];
cx q14[0],q14[1];
y q14[0];
s q14[0];
sdg q14[1];
h q14[1];
x q14[1];
h q14[1];
cx q14[0],q14[1];
sdg q14[0];
h q14[0];
h q14[1];
s q14[1];
cx q14[0],q14[1];
z q14[0];
s q14[0];
h q14[0];
s q14[0];
z q14[1];
cx q14[0],q14[1];
h q14[1];
s q14[1];
y q14[1];
s q14[1];
h q14[1];
s q14[1];
cx q14[0],q14[1];
sdg q14[0];
h q14[0];
sdg q14[1];
h q14[1];
cx q14[0],q14[1];
h q14[0];
s q14[0];
z q14[0];
h q14[0];
h q14[1];
s q14[1];
cx q14[0],q14[1];
h q14[0];
s q14[0];
h q14[0];
s q14[0];
sdg q14[1];
h q14[1];
x q14[1];
s q14[1];
cx q14[0],q14[1];
h q14[0];
s q14[0];
y q14[0];
h q14[0];
s q14[0];
h q14[1];
s q14[1];
y q14[1];
s q14[1];
h q14[1];
s q14[1];
cx q14[0],q14[1];
sdg q14[0];
h q14[0];
sdg q14[1];
h q14[1];
cx q14[0],q14[1];
y q14[0];
s q14[0];
h q14[0];
s q14[0];
y q14[1];
s q14[1];
h q14[1];
s q14[1];
cx q14[0],q14[1];
sdg q14[0];
h q14[0];
y q14[0];
h q14[0];
x q14[1];
s q14[1];
h q14[1];
s q14[1];
cx q14[0],q14[1];
sdg q14[0];
h q14[0];
x q14[0];
s q14[0];
sdg q14[1];
h q14[1];
y q14[1];
cx q14[0],q14[1];
sdg q14[0];
h q14[0];
sdg q14[1];
h q14[1];
cx q14[0],q14[1];
x q14[0];
s q14[0];
z q14[1];
s q14[1];
h q14[1];
s q14[1];
cx q14[0],q14[1];
sdg q14[0];
h q14[0];
sdg q14[1];
h q14[1];
cx q14[0],q14[1];
h q14[0];
s q14[0];
y q14[0];
sdg q14[0];
h q14[0];
z q14[1];
cx q14[0],q14[1];
h q14[0];
s q14[0];
y q14[0];
s q14[0];
h q14[1];
s q14[1];
y q14[1];
s q14[1];
h q14[1];
s q14[1];
cx q14[0],q14[1];
sdg q14[0];
h q14[0];
sdg q14[1];
h q14[1];
cx q14[0],q14[1];
y q14[0];
s q14[0];
h q14[0];
s q14[0];
h q14[1];
s q14[1];
y q14[1];
h q14[1];
cx q14[0],q14[1];
sdg q14[0];
h q14[0];
sdg q14[1];
h q14[1];
cx q14[0],q14[1];
h q14[0];
s q14[0];
h q14[1];
s q14[1];
cx q14[0],q14[1];
x q14[0];
x q14[1];
cx q14[0],q14[1];
h q14[0];
s q14[0];
h q14[0];
h q14[1];
s q14[1];
y q14[1];
h q14[1];
cx q14[0],q14[1];
sdg q14[0];
h q14[0];
y q14[0];
h q14[0];
s q14[0];
h q14[1];
s q14[1];
x q14[1];
s q14[1];
h q14[1];
s q14[1];
cx q14[0],q14[1];
sdg q14[0];
h q14[0];
sdg q14[1];
h q14[1];
cx q14[0],q14[1];
h q14[0];
s q14[0];
h q14[1];
s q14[1];
cx q14[0],q14[1];
z q14[0];
s q14[0];
h q14[0];
s q14[0];
y q14[1];
h q14[1];
s q14[1];
cx q14[0],q14[1];
sdg q14[0];
h q14[0];
y q14[0];
h q14[0];
s q14[0];
z q14[1];
h q14[1];
s q14[1];
cx q14[0],q14[1];
h q14[0];
s q14[0];
h q14[1];
s q14[1];
cx q14[0],q14[1];
s q14[0];
h q14[0];
s q14[0];
x q14[1];
cx q14[0],q14[1];
h q14[0];
s q14[0];
y q14[0];
sdg q14[1];
h q14[1];
sdg q14[1];
h q14[1];
cx q14[0],q14[1];
sdg q14[0];
h q14[0];
sdg q14[0];
h q14[0];
h q14[1];
s q14[1];
cx q14[0],q14[1];
sdg q14[0];
h q14[0];
sdg q14[1];
h q14[1];
cx q14[0],q14[1];
h q14[0];
s q14[0];
h q14[1];
s q14[1];
cx q14[0],q14[1];
x q14[0];
s q14[0];
h q14[0];
s q14[0];
x q14[1];
s q14[1];
h q14[1];
s q14[1];
cx q14[0],q14[1];
sdg q14[0];
h q14[0];
z q14[0];
sdg q14[0];
h q14[0];
sdg q14[1];
h q14[1];
z q14[1];
h q14[1];
s q14[1];
cx q14[0],q14[1];
h q14[0];
s q14[0];
sdg q14[1];
h q14[1];
cx q14[0],q14[1];
y q14[0];
x q14[1];
cx q14[0],q14[1];
y q14[0];
s q14[0];
h q14[1];
s q14[1];
h q14[1];
s q14[1];
cx q14[0],q14[1];
h q14[0];
s q14[0];
sdg q14[1];
h q14[1];
cx q14[0],q14[1];
y q14[0];
h q14[1];
s q14[1];
x q14[1];
s q14[1];
cx q14[0],q14[1];
h q14[0];
s q14[0];
sdg q14[1];
h q14[1];
cx q14[0],q14[1];
z q14[0];
z q14[1];
h q14[1];
s q14[1];
cx q14[0],q14[1];
sdg q14[0];
h q14[0];
y q14[0];
s q14[0];
h q14[0];
s q14[0];
h q14[1];
s q14[1];
z q14[1];
sdg q14[1];
h q14[1];
cx q14[0],q14[1];
sdg q14[0];
h q14[0];
sdg q14[1];
h q14[1];
cx q14[0],q14[1];
h q14[0];
s q14[0];
h q14[1];
s q14[1];
cx q14[0],q14[1];
z q14[0];
h q14[0];
h q14[0];
s q14[0];
z q14[1];
s q14[1];
z q14[1];
s q14[1];
h q14[1];
s q14[1];
cx q14[0],q14[1];
h q14[0];
s q14[0];
z q14[0];
s q14[0];
h q14[0];
s q14[0];
y q14[1];
sdg q14[1];
h q14[1];
cx q14[0],q14[1];
sdg q14[0];
h q14[0];
x q14[0];
y q14[0];
h q14[0];
s q14[0];
h q14[1];
s q14[1];
s q14[1];
z q14[1];
h q14[1];
cx q14[0],q14[1];
s q14[0];
sdg q14[1];
h q14[1];
y q14[1];
h q14[1];
cx q14[0],q14[1];
sdg q14[0];
h q14[0];
h q14[0];
s q14[0];
sdg q14[1];
h q14[1];
y q14[1];
cx q14[0],q14[1];
sdg q14[0];
h q14[0];
x q14[0];
h q14[0];
x q14[1];
cx q14[0],q14[1];
h q14[0];
s q14[0];
sdg q14[1];
h q14[1];
cx q14[0],q14[1];
y q14[0];
s q14[0];
h q14[0];
s q14[0];
h q14[1];
s q14[1];
z q14[1];
h q14[1];
s q14[1];
cx q14[0],q14[1];
h q14[0];
s q14[0];
sdg q14[0];
h q14[0];
sdg q14[1];
h q14[1];
h q14[1];
cx q14[0],q14[1];
h q14[0];
s q14[0];
z q14[0];
sdg q14[0];
h q14[0];
sdg q14[1];
h q14[1];
x q14[1];
s q14[1];
h q14[1];
s q14[1];
cx q14[0],q14[1];
sdg q14[0];
h q14[0];
h q14[1];
s q14[1];
cx q14[0],q14[1];
x q14[0];
sdg q14[0];
h q14[0];
z q14[1];
sdg q14[1];
h q14[1];
cx q14[0],q14[1];
sdg q14[0];
h q14[0];
sdg q14[1];
h q14[1];
cx q14[0],q14[1];
h q14[0];
s q14[0];
x q14[0];
h q14[1];
s q14[1];
y q14[1];
s q14[1];
h q14[1];
s q14[1];
cx q14[0],q14[1];
sdg q14[0];
h q14[0];
sdg q14[1];
h q14[1];
cx q14[0],q14[1];
h q14[0];
s q14[0];
z q14[0];
sdg q14[0];
h q14[0];
y q14[1];
sdg q14[1];
h q14[1];
cx q14[0],q14[1];
h q14[0];
s q14[0];
x q14[0];
h q14[0];
s q14[0];
sdg q14[1];
h q14[1];
y q14[1];
h q14[1];
s q14[1];
cx q14[0],q14[1];
h q14[0];
s q14[0];
h q14[1];
s q14[1];
cx q14[0],q14[1];
z q14[0];
h q14[0];
s q14[0];
y q14[1];
s q14[1];
h q14[1];
s q14[1];
cx q14[0],q14[1];
sdg q14[0];
h q14[0];
h q14[1];
s q14[1];
cx q14[0],q14[1];
h q14[0];
s q14[0];
h q14[0];
x q14[1];
s q14[1];
cx q14[0],q14[1];
s q14[0];
h q14[0];
s q14[0];
sdg q14[1];
h q14[1];
x q14[1];
h q14[1];
s q14[1];
cx q14[0],q14[1];
sdg q14[0];
h q14[0];
z q14[0];
sdg q14[0];
h q14[0];
sdg q14[1];
h q14[1];
h q14[1];
s q14[1];
cx q14[0],q14[1];
sdg q14[0];
h q14[0];
sdg q14[1];
h q14[1];
cx q14[0],q14[1];
h q14[0];
s q14[0];
x q14[0];
h q14[0];
h q14[1];
s q14[1];
y q14[1];
h q14[1];
s q14[1];
cx q14[0],q14[1];
x q14[0];
s q14[0];
h q14[0];
s q14[0];
y q14[0];
x q14[1];
h q14[1];
s q14[1];
z q14[1];
h q14[1];
s q14[1];
cx q14[0],q14[1];
sdg q14[0];
h q14[0];
h q14[0];
h q14[1];
s q14[1];
y q14[1];
cx q14[0],q14[1];
h q14[0];
s q14[0];
sdg q14[0];
h q14[0];
sdg q14[1];
h q14[1];
z q14[1];
cx q14[0],q14[1];
sdg q14[0];
h q14[0];
sdg q14[1];
h q14[1];
cx q14[0],q14[1];
x q14[0];
h q14[1];
s q14[1];
x q14[1];
sdg q14[1];
h q14[1];
cx q14[0],q14[1];
sdg q14[0];
h q14[0];
sdg q14[1];
h q14[1];
cx q14[0],q14[1];
h q14[0];
s q14[0];
h q14[1];
s q14[1];
cx q14[0],q14[1];
x q14[0];
s q14[0];
x q14[1];
h q14[1];
cx q14[0],q14[1];
sdg q14[0];
h q14[0];
z q14[0];
h q14[0];
h q14[1];
s q14[1];
y q14[1];
sdg q14[1];
h q14[1];
cx q14[0],q14[1];
z q14[0];
x q14[1];
cx q14[0],q14[1];
x q14[0];
sdg q14[0];
h q14[0];
sdg q14[1];
h q14[1];
y q14[1];
cx q14[0],q14[1];
sdg q14[0];
h q14[0];
h q14[1];
s q14[1];
cx q14[0],q14[1];
y q14[0];
h q14[0];
sdg q14[1];
h q14[1];
cx q14[0],q14[1];
h q14[0];
s q14[0];
sdg q14[1];
h q14[1];
cx q14[0],q14[1];
x q14[0];
sdg q14[0];
h q14[0];
h q14[1];
s q14[1];
z q14[1];
sdg q14[1];
h q14[1];
cx q14[0],q14[1];
sdg q14[0];
h q14[0];
y q14[0];
h q14[0];
s q14[0];
h q14[1];
s q14[1];
z q14[1];
s q14[1];
cx q14[0],q14[1];
sdg q14[0];
h q14[0];
h q14[1];
s q14[1];
cx q14[0],q14[1];
h q14[0];
s q14[0];
z q14[0];
s q14[0];
h q14[0];
s q14[0];
y q14[1];
s q14[1];
cx q14[0],q14[1];
h q14[0];
s q14[0];
z q14[0];
sdg q14[0];
h q14[0];
sdg q14[1];
h q14[1];
x q14[1];
cx q14[0],q14[1];
sdg q14[0];
h q14[0];
h q14[1];
s q14[1];
cx q14[0],q14[1];
h q14[0];
s q14[0];
sdg q14[0];
h q14[0];
h q14[1];
cx q14[0],q14[1];
h q14[0];
s q14[0];
x q14[0];
s q14[0];
h q14[0];
s q14[0];
sdg q14[1];
h q14[1];
z q14[1];
sdg q14[1];
h q14[1];
cx q14[0],q14[1];
h q14[0];
s q14[0];
z q14[0];
s q14[0];
sdg q14[1];
h q14[1];
x q14[1];
h q14[1];
cx q14[0],q14[1];
sdg q14[0];
h q14[0];
sdg q14[1];
h q14[1];
cx q14[0],q14[1];
h q14[0];
s q14[0];
z q14[0];
h q14[0];
s q14[0];
z q14[1];
s q14[1];
cx q14[0],q14[1];
sdg q14[0];
h q14[0];
h q14[1];
s q14[1];
cx q14[0],q14[1];
h q14[0];
y q14[1];
h q14[1];
cx q14[0],q14[1];
sdg q14[0];
h q14[0];
z q14[0];
h q14[0];
s q14[0];
sdg q14[1];
h q14[1];
y q14[1];
h q14[1];
s q14[1];
cx q14[0],q14[1];
x q14[0];
h q14[0];
h q14[1];
s q14[1];
x q14[1];
cx q14[0],q14[1];
sdg q14[0];
h q14[0];
h q14[0];
s q14[0];
sdg q14[1];
h q14[1];
h q14[1];
s q14[1];
cx q14[0],q14[1];
sdg q14[1];
h q14[1];
y q14[1];
sdg q14[1];
h q14[1];
cx q14[0],q14[1];
h q14[0];
s q14[0];
h q14[1];
s q14[1];
cx q14[0],q14[1];
h q14[0];
x q14[1];
s q14[1];
cx q14[0],q14[1];
h q14[0];
s q14[0];
sdg q14[1];
h q14[1];
cx q14[0],q14[1];
x q14[0];
s q14[0];
h q14[0];
s q14[0];
h q14[1];
s q14[1];
z q14[1];
s q14[1];
cx q14[0],q14[1];
x q14[0];
s q14[0];
h q14[0];
s q14[0];
h q14[1];
s q14[1];
z q14[1];
s q14[1];
cx q14[0],q14[1];
h q14[0];
s q14[0];
z q14[0];
s q14[0];
h q14[0];
s q14[0];
sdg q14[1];
h q14[1];
y q14[1];
cx q14[0],q14[1];
sdg q14[0];
h q14[0];
x q14[0];
s q14[0];
h q14[0];
s q14[0];
sdg q14[1];
h q14[1];
y q14[1];
h q14[1];
s q14[1];
cx q14[0],q14[1];
sdg q14[0];
h q14[0];
h q14[1];
s q14[1];
cx q14[0],q14[1];
h q14[0];
s q14[0];
y q14[0];
s q14[0];
h q14[1];
cx q14[0],q14[1];
h q14[0];
s q14[0];
sdg q14[1];
h q14[1];
cx q14[0],q14[1];
z q14[0];
h q14[0];
z q14[1];
cx q14[0],q14[1];
h q14[0];
s q14[0];
h q14[1];
s q14[1];
cx q14[0],q14[1];
x q14[0];
y q14[1];
s q14[1];
h q14[1];
s q14[1];
cx q14[0],q14[1];
sdg q14[0];
h q14[0];
s q14[0];
sdg q14[1];
h q14[1];
cx q14[0],q14[1];
sdg q14[0];
h q14[0];
sdg q14[1];
h q14[1];
cx q14[0],q14[1];
h q14[0];
s q14[0];
h q14[1];
s q14[1];
cx q14[0],q14[1];
z q14[0];
s q14[0];
h q14[0];
s q14[0];
y q14[1];
s q14[1];
cx q14[0],q14[1];
sdg q14[0];
h q14[0];
sdg q14[1];
h q14[1];
cx q14[0],q14[1];
h q14[0];
s q14[0];
x q14[0];
h q14[1];
s q14[1];
s q14[1];
cx q14[0],q14[1];
sdg q14[0];
h q14[0];
h q14[1];
s q14[1];
cx q14[0],q14[1];
y q14[0];
h q14[0];
s q14[0];
h q14[1];
s q14[1];
cx q14[0],q14[1];
h q14[0];
s q14[0];
s q14[0];
x q14[0];
h q14[0];
s q14[0];
h q14[1];
s q14[1];
z q14[1];
s q14[1];
h q14[1];
s q14[1];
s q14[1];
cx q14[0],q14[1];
h q14[0];
s q14[0];
x q14[0];
h q14[1];
s q14[1];
z q14[1];
cx q14[0],q14[1];
h q14[0];
s q14[0];
sdg q14[1];
h q14[1];
cx q14[0],q14[1];
s q14[0];
h q14[0];
s q14[0];
x q14[1];
s q14[1];
cx q14[0],q14[1];
y q14[0];
h q14[0];
h q14[1];
s q14[1];
sdg q14[1];
h q14[1];
cx q14[0],q14[1];
sdg q14[0];
h q14[0];
sdg q14[1];
h q14[1];
cx q14[0],q14[1];
h q14[0];
s q14[0];
h q14[1];
s q14[1];
cx q14[0],q14[1];
x q14[1];
sdg q14[1];
h q14[1];
cx q14[0],q14[1];
sdg q14[0];
h q14[0];
z q14[0];
h q14[0];
h q14[1];
s q14[1];
x q14[1];
sdg q14[1];
h q14[1];
cx q14[0],q14[1];
sdg q14[0];
h q14[0];
sdg q14[1];
h q14[1];
cx q14[0],q14[1];
x q14[0];
s q14[0];
x q14[1];
h q14[1];
cx q14[0],q14[1];
x q14[0];
s q14[0];
h q14[0];
s q14[0];
x q14[1];
sdg q14[1];
h q14[1];
cx q14[0],q14[1];
sdg q14[0];
h q14[0];
z q14[0];
h q14[0];
s q14[0];
sdg q14[1];
h q14[1];
z q14[1];
s q14[1];
h q14[1];
s q14[1];
cx q14[0],q14[1];
z q14[0];
sdg q14[1];
h q14[1];
x q14[1];
h q14[1];
s q14[1];
cx q14[0],q14[1];
y q14[0];
sdg q14[0];
h q14[0];
y q14[1];
s q14[1];
cx q14[0],q14[1];
sdg q14[0];
h q14[0];
x q14[0];
s q14[0];
y q14[1];
s q14[1];
cx q14[0],q14[1];
sdg q14[0];
h q14[0];
sdg q14[1];
h q14[1];
cx q14[0],q14[1];
h q14[0];
s q14[0];
z q14[0];
h q14[0];
x q14[1];
s q14[1];
cx q14[0],q14[1];
sdg q14[0];
h q14[0];
h q14[1];
s q14[1];
cx q14[0],q14[1];
h q14[0];
s q14[0];
x q14[1];
sdg q14[1];
h q14[1];
cx q14[0],q14[1];
sdg q14[0];
h q14[0];
sdg q14[1];
h q14[1];
cx q14[0],q14[1];
h q14[0];
s q14[0];
h q14[1];
s q14[1];
cx q14[0],q14[1];
x q14[0];
h q14[0];
s q14[0];
x q14[1];
h q14[1];
cx q14[0],q14[1];
sdg q14[0];
h q14[0];
sdg q14[1];
h q14[1];
cx q14[0],q14[1];
h q14[0];
s q14[0];
x q14[0];
sdg q14[0];
h q14[0];
y q14[1];
h q14[1];
s q14[1];
cx q14[0],q14[1];
h q14[0];
s q14[0];
sdg q14[1];
h q14[1];
