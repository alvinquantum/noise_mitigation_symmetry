OPENQASM 2.0;
include "qelib1.inc";
qreg q24[2];
z q24[0];
sdg q24[0];
h q24[0];
z q24[1];
h q24[1];
s q24[1];
cx q24[0],q24[1];
sdg q24[0];
h q24[0];
sdg q24[1];
h q24[1];
cx q24[0],q24[1];
h q24[0];
s q24[0];
h q24[1];
s q24[1];
cx q24[0],q24[1];
h q24[1];
cx q24[0],q24[1];
sdg q24[0];
h q24[0];
s q24[0];
h q24[1];
s q24[1];
y q24[1];
h q24[1];
s q24[1];
cx q24[0],q24[1];
sdg q24[0];
h q24[0];
sdg q24[1];
h q24[1];
cx q24[0],q24[1];
x q24[0];
s q24[0];
y q24[0];
h q24[0];
s q24[0];
h q24[1];
s q24[1];
y q24[1];
s q24[1];
s q24[1];
cx q24[0],q24[1];
s q24[0];
h q24[0];
s q24[0];
h q24[1];
cx q24[0],q24[1];
sdg q24[0];
h q24[0];
h q24[1];
s q24[1];
cx q24[0],q24[1];
y q24[0];
sdg q24[0];
h q24[0];
x q24[1];
h q24[1];
s q24[1];
cx q24[0],q24[1];
sdg q24[0];
h q24[0];
sdg q24[0];
h q24[0];
sdg q24[1];
h q24[1];
y q24[1];
cx q24[0],q24[1];
y q24[0];
h q24[0];
s q24[0];
x q24[1];
s q24[1];
cx q24[0],q24[1];
sdg q24[0];
h q24[0];
sdg q24[1];
h q24[1];
cx q24[0],q24[1];
h q24[0];
s q24[0];
s q24[0];
h q24[1];
s q24[1];
z q24[1];
h q24[1];
s q24[1];
cx q24[0],q24[1];
z q24[0];
h q24[0];
s q24[0];
y q24[1];
h q24[1];
cx q24[0],q24[1];
sdg q24[0];
h q24[0];
sdg q24[1];
h q24[1];
cx q24[0],q24[1];
h q24[0];
s q24[0];
h q24[1];
s q24[1];
cx q24[0],q24[1];
y q24[0];
s q24[0];
h q24[0];
s q24[0];
y q24[1];
h q24[1];
s q24[1];
cx q24[0],q24[1];
y q24[0];
sdg q24[0];
h q24[0];
sdg q24[1];
h q24[1];
y q24[1];
h q24[1];
s q24[1];
cx q24[0],q24[1];
z q24[0];
h q24[1];
s q24[1];
h q24[1];
cx q24[0],q24[1];
h q24[0];
s q24[0];
s q24[0];
sdg q24[1];
h q24[1];
x q24[1];
cx q24[0],q24[1];
sdg q24[0];
h q24[0];
sdg q24[1];
h q24[1];
cx q24[0],q24[1];
h q24[0];
s q24[0];
z q24[0];
s q24[0];
h q24[1];
cx q24[0],q24[1];
h q24[0];
s q24[0];
z q24[0];
h q24[0];
s q24[0];
h q24[1];
s q24[1];
s q24[1];
h q24[1];
s q24[1];
cx q24[0],q24[1];
sdg q24[0];
h q24[0];
sdg q24[1];
h q24[1];
cx q24[0],q24[1];
h q24[0];
s q24[0];
y q24[0];
h q24[0];
cx q24[0],q24[1];
sdg q24[0];
h q24[0];
sdg q24[1];
h q24[1];
cx q24[0],q24[1];
h q24[0];
s q24[0];
h q24[1];
s q24[1];
y q24[1];
h q24[1];
cx q24[0],q24[1];
sdg q24[0];
h q24[0];
sdg q24[1];
h q24[1];
cx q24[0],q24[1];
h q24[0];
s q24[0];
x q24[0];
s q24[0];
h q24[0];
s q24[0];
h q24[1];
s q24[1];
s q24[1];
h q24[1];
s q24[1];
cx q24[0],q24[1];
sdg q24[0];
h q24[0];
h q24[1];
s q24[1];
cx q24[0],q24[1];
x q24[0];
h q24[0];
s q24[0];
z q24[1];
cx q24[0],q24[1];
sdg q24[0];
h q24[0];
sdg q24[1];
h q24[1];
cx q24[0],q24[1];
h q24[0];
s q24[0];
z q24[0];
h q24[0];
s q24[0];
x q24[1];
cx q24[0],q24[1];
sdg q24[0];
h q24[0];
sdg q24[1];
h q24[1];
cx q24[0],q24[1];
s q24[0];
h q24[0];
s q24[0];
h q24[1];
s q24[1];
x q24[1];
sdg q24[1];
h q24[1];
cx q24[0],q24[1];
sdg q24[0];
h q24[0];
sdg q24[1];
h q24[1];
cx q24[0],q24[1];
h q24[0];
s q24[0];
h q24[1];
s q24[1];
cx q24[0],q24[1];
x q24[0];
h q24[0];
h q24[1];
s q24[1];
cx q24[0],q24[1];
sdg q24[0];
h q24[0];
sdg q24[1];
h q24[1];
cx q24[0],q24[1];
h q24[0];
s q24[0];
h q24[1];
s q24[1];
cx q24[0],q24[1];
z q24[0];
h q24[0];
y q24[1];
h q24[1];
cx q24[0],q24[1];
h q24[0];
s q24[0];
h q24[0];
s q24[0];
sdg q24[1];
h q24[1];
x q24[1];
s q24[1];
h q24[1];
s q24[1];
cx q24[0],q24[1];
sdg q24[0];
h q24[0];
sdg q24[1];
h q24[1];
cx q24[0],q24[1];
h q24[0];
s q24[0];
h q24[1];
s q24[1];
cx q24[0],q24[1];
y q24[0];
s q24[0];
y q24[1];
cx q24[0],q24[1];
sdg q24[0];
h q24[0];
z q24[0];
s q24[0];
h q24[0];
s q24[0];
sdg q24[1];
h q24[1];
z q24[1];
h q24[1];
cx q24[0],q24[1];
h q24[0];
s q24[0];
x q24[0];
h q24[0];
s q24[0];
s q24[1];
cx q24[0],q24[1];
sdg q24[0];
h q24[0];
y q24[0];
h q24[0];
s q24[0];
sdg q24[1];
h q24[1];
x q24[1];
h q24[1];
cx q24[0],q24[1];
sdg q24[0];
h q24[0];
sdg q24[1];
h q24[1];
cx q24[0],q24[1];
h q24[0];
s q24[0];
z q24[0];
x q24[1];
s q24[1];
h q24[1];
s q24[1];
cx q24[0],q24[1];
sdg q24[0];
h q24[0];
sdg q24[1];
h q24[1];
cx q24[0],q24[1];
x q24[0];
h q24[0];
s q24[0];
x q24[0];
h q24[1];
s q24[1];
s q24[1];
h q24[1];
s q24[1];
h q24[1];
s q24[1];
cx q24[0],q24[1];
sdg q24[0];
h q24[0];
h q24[0];
y q24[1];
s q24[1];
h q24[1];
s q24[1];
cx q24[0],q24[1];
sdg q24[0];
h q24[0];
h q24[1];
s q24[1];
cx q24[0],q24[1];
z q24[0];
sdg q24[0];
h q24[0];
z q24[1];
h q24[1];
cx q24[0],q24[1];
z q24[0];
y q24[1];
sdg q24[1];
h q24[1];
cx q24[0],q24[1];
sdg q24[0];
h q24[0];
z q24[0];
h q24[0];
s q24[0];
x q24[0];
sdg q24[0];
h q24[0];
h q24[1];
s q24[1];
y q24[1];
s q24[1];
h q24[1];
s q24[1];
sdg q24[1];
h q24[1];
cx q24[0],q24[1];
sdg q24[0];
h q24[0];
sdg q24[0];
h q24[0];
z q24[0];
s q24[0];
sdg q24[1];
h q24[1];
y q24[1];
h q24[1];
cx q24[0],q24[1];
sdg q24[0];
h q24[0];
h q24[1];
s q24[1];
cx q24[0],q24[1];
y q24[0];
h q24[0];
s q24[0];
x q24[1];
h q24[1];
cx q24[0],q24[1];
x q24[0];
s q24[0];
sdg q24[1];
h q24[1];
z q24[1];
h q24[1];
s q24[1];
cx q24[0],q24[1];
h q24[0];
s q24[0];
sdg q24[1];
h q24[1];
cx q24[0],q24[1];
x q24[0];
sdg q24[0];
h q24[0];
h q24[1];
s q24[1];
z q24[1];
s q24[1];
h q24[1];
s q24[1];
cx q24[0],q24[1];
sdg q24[0];
h q24[0];
z q24[0];
sdg q24[0];
h q24[0];
sdg q24[1];
h q24[1];
h q24[1];
cx q24[0],q24[1];
sdg q24[0];
h q24[0];
sdg q24[1];
h q24[1];
cx q24[0],q24[1];
h q24[0];
s q24[0];
h q24[1];
s q24[1];
cx q24[0],q24[1];
y q24[0];
s q24[0];
z q24[1];
sdg q24[1];
h q24[1];
cx q24[0],q24[1];
z q24[0];
h q24[0];
s q24[0];
z q24[0];
s q24[0];
x q24[1];
sdg q24[1];
h q24[1];
y q24[1];
sdg q24[1];
h q24[1];
cx q24[0],q24[1];
y q24[0];
s q24[0];
y q24[1];
sdg q24[1];
h q24[1];
cx q24[0],q24[1];
sdg q24[0];
h q24[0];
sdg q24[1];
h q24[1];
cx q24[0],q24[1];
h q24[0];
s q24[0];
x q24[0];
s q24[0];
h q24[1];
s q24[1];
z q24[1];
h q24[1];
cx q24[0],q24[1];
h q24[0];
s q24[0];
z q24[0];
s q24[0];
z q24[1];
sdg q24[1];
h q24[1];
cx q24[0],q24[1];
sdg q24[0];
h q24[0];
sdg q24[1];
h q24[1];
cx q24[0],q24[1];
h q24[0];
s q24[0];
h q24[1];
s q24[1];
cx q24[0],q24[1];
x q24[0];
sdg q24[0];
h q24[0];
cx q24[0],q24[1];
sdg q24[0];
h q24[0];
sdg q24[1];
h q24[1];
cx q24[0],q24[1];
x q24[0];
h q24[0];
cx q24[0],q24[1];
sdg q24[0];
h q24[0];
x q24[0];
s q24[0];
h q24[0];
s q24[0];
y q24[1];
cx q24[0],q24[1];
sdg q24[0];
h q24[0];
s q24[0];
h q24[1];
s q24[1];
h q24[1];
s q24[1];
cx q24[0],q24[1];
z q24[0];
sdg q24[0];
h q24[0];
sdg q24[1];
h q24[1];
h q24[1];
cx q24[0],q24[1];
h q24[0];
s q24[0];
z q24[0];
s q24[0];
sdg q24[1];
h q24[1];
x q24[1];
s q24[1];
h q24[1];
s q24[1];
cx q24[0],q24[1];
sdg q24[0];
h q24[0];
y q24[0];
h q24[0];
s q24[0];
h q24[1];
s q24[1];
x q24[1];
cx q24[0],q24[1];
sdg q24[0];
h q24[0];
sdg q24[1];
h q24[1];
cx q24[0],q24[1];
h q24[0];
s q24[0];
z q24[0];
h q24[0];
z q24[1];
s q24[1];
h q24[1];
s q24[1];
cx q24[0],q24[1];
sdg q24[0];
h q24[0];
x q24[0];
sdg q24[1];
h q24[1];
h q24[1];
s q24[1];
cx q24[0],q24[1];
h q24[0];
s q24[0];
y q24[0];
cx q24[0],q24[1];
sdg q24[0];
h q24[0];
sdg q24[1];
h q24[1];
cx q24[0],q24[1];
y q24[0];
h q24[0];
s q24[0];
x q24[1];
h q24[1];
sdg q24[1];
h q24[1];
cx q24[0],q24[1];
y q24[0];
h q24[1];
s q24[1];
x q24[1];
h q24[1];
s q24[1];
cx q24[0],q24[1];
sdg q24[0];
h q24[0];
sdg q24[1];
h q24[1];
cx q24[0],q24[1];
h q24[0];
s q24[0];
h q24[1];
s q24[1];
cx q24[0],q24[1];
y q24[0];
h q24[0];
s q24[0];
y q24[1];
s q24[1];
h q24[1];
s q24[1];
cx q24[0],q24[1];
sdg q24[0];
h q24[0];
sdg q24[1];
h q24[1];
cx q24[0],q24[1];
h q24[0];
s q24[0];
z q24[0];
s q24[0];
h q24[0];
s q24[0];
sdg q24[1];
h q24[1];
cx q24[0],q24[1];
sdg q24[0];
h q24[0];
sdg q24[1];
h q24[1];
cx q24[0],q24[1];
h q24[0];
s q24[0];
y q24[0];
s q24[0];
h q24[0];
s q24[0];
x q24[1];
h q24[1];
cx q24[0],q24[1];
sdg q24[0];
h q24[0];
x q24[0];
s q24[0];
h q24[1];
s q24[1];
z q24[1];
sdg q24[1];
h q24[1];
cx q24[0],q24[1];
h q24[0];
s q24[0];
h q24[1];
s q24[1];
cx q24[0],q24[1];
s q24[0];
z q24[1];
sdg q24[1];
h q24[1];
cx q24[0],q24[1];
sdg q24[0];
h q24[0];
sdg q24[1];
h q24[1];
cx q24[0],q24[1];
h q24[0];
s q24[0];
h q24[1];
s q24[1];
cx q24[0],q24[1];
sdg q24[0];
h q24[0];
y q24[1];
s q24[1];
cx q24[0],q24[1];
sdg q24[0];
h q24[0];
x q24[0];
sdg q24[1];
h q24[1];
y q24[1];
s q24[1];
h q24[1];
s q24[1];
cx q24[0],q24[1];
sdg q24[0];
h q24[0];
h q24[1];
s q24[1];
cx q24[0],q24[1];
z q24[0];
h q24[1];
s q24[1];
cx q24[0],q24[1];
sdg q24[0];
h q24[0];
sdg q24[1];
h q24[1];
cx q24[0],q24[1];
x q24[0];
h q24[1];
cx q24[0],q24[1];
sdg q24[0];
h q24[0];
y q24[0];
h q24[0];
sdg q24[1];
h q24[1];
z q24[1];
h q24[1];
cx q24[0],q24[1];
sdg q24[0];
h q24[0];
y q24[0];
h q24[0];
s q24[0];
h q24[1];
s q24[1];
y q24[1];
cx q24[0],q24[1];
sdg q24[0];
h q24[0];
sdg q24[1];
h q24[1];
cx q24[0],q24[1];
h q24[0];
s q24[0];
x q24[0];
h q24[0];
h q24[1];
s q24[1];
s q24[1];
cx q24[0],q24[1];
sdg q24[0];
h q24[0];
s q24[0];
h q24[0];
s q24[0];
h q24[1];
s q24[1];
h q24[1];
s q24[1];
cx q24[0],q24[1];
sdg q24[0];
h q24[0];
h q24[1];
s q24[1];
cx q24[0],q24[1];
h q24[0];
s q24[0];
x q24[0];
h q24[0];
s q24[1];
cx q24[0],q24[1];
sdg q24[0];
h q24[0];
sdg q24[1];
h q24[1];
cx q24[0],q24[1];
h q24[0];
s q24[0];
x q24[0];
s q24[0];
h q24[1];
s q24[1];
y q24[1];
s q24[1];
cx q24[0],q24[1];
sdg q24[0];
h q24[0];
sdg q24[1];
h q24[1];
cx q24[0],q24[1];
x q24[0];
s q24[0];
h q24[1];
s q24[1];
z q24[1];
cx q24[0],q24[1];
sdg q24[0];
h q24[0];
sdg q24[1];
h q24[1];
cx q24[0],q24[1];
z q24[0];
h q24[0];
h q24[1];
s q24[1];
s q24[1];
cx q24[0],q24[1];
h q24[0];
s q24[0];
sdg q24[1];
h q24[1];
cx q24[0],q24[1];
y q24[0];
sdg q24[0];
h q24[0];
x q24[1];
cx q24[0],q24[1];
h q24[0];
s q24[0];
sdg q24[1];
h q24[1];
cx q24[0],q24[1];
z q24[0];
h q24[0];
s q24[0];
h q24[1];
s q24[1];
z q24[1];
h q24[1];
cx q24[0],q24[1];
h q24[0];
s q24[0];
h q24[1];
s q24[1];
cx q24[0],q24[1];
s q24[0];
h q24[0];
s q24[0];
s q24[1];
h q24[1];
s q24[1];
cx q24[0],q24[1];
h q24[0];
s q24[0];
z q24[0];
h q24[1];
s q24[1];
y q24[1];
s q24[1];
cx q24[0],q24[1];
h q24[1];
s q24[1];
y q24[1];
h q24[1];
cx q24[0],q24[1];
sdg q24[0];
h q24[0];
sdg q24[1];
h q24[1];
cx q24[0],q24[1];
h q24[0];
s q24[0];
h q24[1];
s q24[1];
z q24[1];
s q24[1];
h q24[1];
s q24[1];
cx q24[0],q24[1];
y q24[0];
h q24[0];
s q24[0];
z q24[1];
h q24[1];
cx q24[0],q24[1];
sdg q24[0];
h q24[0];
x q24[0];
s q24[0];
h q24[0];
s q24[0];
y q24[1];
h q24[1];
cx q24[0],q24[1];
z q24[0];
x q24[1];
cx q24[0],q24[1];
sdg q24[0];
h q24[0];
sdg q24[1];
h q24[1];
cx q24[0],q24[1];
h q24[0];
s q24[0];
h q24[0];
s q24[0];
x q24[1];
s q24[1];
h q24[1];
s q24[1];
cx q24[0],q24[1];
s q24[0];
sdg q24[1];
h q24[1];
y q24[1];
h q24[1];
s q24[1];
cx q24[0],q24[1];
h q24[0];
s q24[0];
x q24[0];
sdg q24[0];
h q24[0];
h q24[1];
s q24[1];
s q24[1];
h q24[1];
s q24[1];
cx q24[0],q24[1];
sdg q24[0];
h q24[0];
sdg q24[1];
h q24[1];
cx q24[0],q24[1];
h q24[0];
s q24[0];
y q24[0];
s q24[0];
x q24[1];
s q24[1];
cx q24[0],q24[1];
sdg q24[0];
h q24[0];
x q24[0];
h q24[0];
s q24[0];
h q24[1];
s q24[1];
x q24[1];
sdg q24[1];
h q24[1];
cx q24[0],q24[1];
h q24[0];
s q24[0];
sdg q24[1];
h q24[1];
cx q24[0],q24[1];
z q24[0];
sdg q24[0];
h q24[0];
x q24[1];
s q24[1];
h q24[1];
s q24[1];
cx q24[0],q24[1];
sdg q24[0];
h q24[0];
sdg q24[1];
h q24[1];
cx q24[0],q24[1];
h q24[0];
s q24[0];
h q24[0];
y q24[1];
s q24[1];
cx q24[0],q24[1];
sdg q24[0];
h q24[0];
x q24[0];
s q24[0];
y q24[0];
sdg q24[0];
h q24[0];
sdg q24[1];
h q24[1];
z q24[1];
h q24[1];
s q24[1];
x q24[1];
s q24[1];
h q24[1];
s q24[1];
cx q24[0],q24[1];
sdg q24[0];
h q24[0];
y q24[0];
sdg q24[1];
h q24[1];
cx q24[0],q24[1];
sdg q24[0];
h q24[0];
h q24[1];
s q24[1];
cx q24[0],q24[1];
h q24[0];
s q24[0];
y q24[0];
z q24[1];
s q24[1];
cx q24[0],q24[1];
sdg q24[0];
h q24[0];
sdg q24[1];
h q24[1];
cx q24[0],q24[1];
x q24[0];
h q24[0];
s q24[0];
sdg q24[1];
h q24[1];
cx q24[0],q24[1];
sdg q24[0];
h q24[0];
h q24[1];
s q24[1];
cx q24[0],q24[1];
h q24[0];
s q24[0];
z q24[0];
s q24[0];
h q24[0];
s q24[0];
s q24[1];
cx q24[0],q24[1];
sdg q24[0];
h q24[0];
h q24[1];
s q24[1];
cx q24[0],q24[1];
x q24[0];
s q24[0];
h q24[0];
s q24[0];
z q24[1];
s q24[1];
h q24[1];
s q24[1];
cx q24[0],q24[1];
sdg q24[0];
h q24[0];
h q24[0];
s q24[1];
cx q24[0],q24[1];
sdg q24[0];
h q24[0];
x q24[0];
h q24[0];
x q24[1];
s q24[1];
h q24[1];
s q24[1];
cx q24[0],q24[1];
sdg q24[0];
h q24[0];
h q24[1];
s q24[1];
cx q24[0],q24[1];
h q24[0];
s q24[0];
h q24[1];
cx q24[0],q24[1];
h q24[0];
s q24[0];
sdg q24[1];
h q24[1];
cx q24[0],q24[1];
s q24[0];
h q24[1];
s q24[1];
x q24[1];
s q24[1];
h q24[1];
s q24[1];
cx q24[0],q24[1];
sdg q24[0];
h q24[0];
sdg q24[1];
h q24[1];
cx q24[0],q24[1];
h q24[0];
s q24[0];
sdg q24[0];
h q24[0];
y q24[1];
sdg q24[1];
h q24[1];
cx q24[0],q24[1];
h q24[0];
s q24[0];
z q24[0];
s q24[0];
sdg q24[1];
h q24[1];
h q24[1];
cx q24[0],q24[1];
sdg q24[0];
h q24[0];
h q24[1];
s q24[1];
cx q24[0],q24[1];
h q24[0];
s q24[0];
z q24[0];
x q24[1];
s q24[1];
cx q24[0],q24[1];
h q24[0];
s q24[0];
z q24[0];
h q24[0];
s q24[0];
sdg q24[1];
h q24[1];
z q24[1];
s q24[1];
cx q24[0],q24[1];
h q24[0];
s q24[0];
h q24[1];
s q24[1];
cx q24[0],q24[1];
x q24[0];
s q24[0];
y q24[1];
h q24[1];
s q24[1];
cx q24[0],q24[1];
sdg q24[0];
h q24[0];
h q24[1];
s q24[1];
cx q24[0],q24[1];
h q24[0];
s q24[0];
x q24[0];
s q24[0];
h q24[0];
s q24[0];
sdg q24[1];
h q24[1];
cx q24[0],q24[1];
h q24[0];
s q24[0];
sdg q24[1];
h q24[1];
cx q24[0],q24[1];
y q24[0];
s q24[0];
h q24[0];
s q24[0];
y q24[1];
sdg q24[1];
h q24[1];
cx q24[0],q24[1];
sdg q24[0];
h q24[0];
h q24[1];
s q24[1];
cx q24[0],q24[1];
sdg q24[0];
h q24[0];
y q24[1];
s q24[1];
h q24[1];
s q24[1];
cx q24[0],q24[1];
h q24[0];
s q24[0];
x q24[0];
s q24[0];
h q24[0];
s q24[0];
z q24[0];
s q24[0];
sdg q24[1];
h q24[1];
h q24[1];
s q24[1];
s q24[1];
h q24[1];
s q24[1];
cx q24[0],q24[1];
sdg q24[0];
h q24[0];
sdg q24[1];
h q24[1];
cx q24[0],q24[1];
h q24[0];
s q24[0];
h q24[1];
s q24[1];
cx q24[0],q24[1];
y q24[0];
h q24[0];
s q24[0];
z q24[1];
sdg q24[1];
h q24[1];
cx q24[0],q24[1];
sdg q24[0];
h q24[0];
h q24[0];
z q24[1];
sdg q24[1];
h q24[1];
cx q24[0],q24[1];
h q24[0];
s q24[0];
y q24[0];
z q24[1];
sdg q24[1];
h q24[1];
cx q24[0],q24[1];
sdg q24[0];
h q24[0];
y q24[0];
sdg q24[0];
h q24[0];
h q24[1];
s q24[1];
sdg q24[1];
h q24[1];
cx q24[0],q24[1];
sdg q24[0];
h q24[0];
z q24[1];
s q24[1];
cx q24[0],q24[1];
sdg q24[0];
h q24[0];
h q24[0];
s q24[0];
s q24[1];
h q24[1];
s q24[1];
cx q24[0],q24[1];
sdg q24[0];
h q24[0];
h q24[1];
s q24[1];
cx q24[0],q24[1];
h q24[0];
s q24[0];
y q24[0];
h q24[0];
cx q24[0],q24[1];
sdg q24[0];
h q24[0];
h q24[1];
s q24[1];
cx q24[0],q24[1];
h q24[0];
s q24[0];
h q24[0];
s q24[0];
s q24[1];
cx q24[0],q24[1];
sdg q24[0];
h q24[0];
sdg q24[1];
h q24[1];
cx q24[0],q24[1];
h q24[0];
s q24[0];
y q24[0];
s q24[0];
h q24[1];
s q24[1];
y q24[1];
cx q24[0],q24[1];
h q24[0];
s q24[0];
h q24[1];
s q24[1];
cx q24[0],q24[1];
z q24[0];
h q24[0];
z q24[1];
cx q24[0],q24[1];
sdg q24[0];
h q24[0];
h q24[1];
s q24[1];
cx q24[0],q24[1];
h q24[0];
s q24[0];
h q24[0];
y q24[1];
sdg q24[1];
h q24[1];
cx q24[0],q24[1];
z q24[0];
s q24[0];
h q24[0];
s q24[0];
x q24[1];
cx q24[0],q24[1];
h q24[0];
s q24[0];
x q24[0];
sdg q24[1];
h q24[1];
h q24[1];
cx q24[0],q24[1];
h q24[0];
s q24[0];
z q24[1];
h q24[1];
s q24[1];
cx q24[0],q24[1];
h q24[0];
s q24[0];
sdg q24[1];
h q24[1];
cx q24[0],q24[1];
s q24[0];
h q24[1];
s q24[1];
z q24[1];
sdg q24[1];
h q24[1];
cx q24[0],q24[1];
h q24[0];
s q24[0];
y q24[0];
s q24[0];
h q24[0];
s q24[0];
h q24[1];
s q24[1];
cx q24[0],q24[1];
sdg q24[0];
h q24[0];
sdg q24[1];
h q24[1];
cx q24[0],q24[1];
h q24[0];
s q24[0];
x q24[0];
s q24[0];
h q24[0];
s q24[0];
s q24[1];
h q24[1];
s q24[1];
cx q24[0],q24[1];
sdg q24[0];
h q24[0];
sdg q24[1];
h q24[1];
cx q24[0],q24[1];
h q24[0];
s q24[0];
y q24[0];
h q24[1];
s q24[1];
z q24[1];
s q24[1];
cx q24[0],q24[1];
h q24[0];
s q24[0];
h q24[1];
s q24[1];
cx q24[0],q24[1];
y q24[0];
s q24[0];
y q24[1];
sdg q24[1];
h q24[1];
cx q24[0],q24[1];
h q24[0];
s q24[0];
x q24[0];
s q24[0];
h q24[0];
s q24[0];
h q24[1];
s q24[1];
z q24[1];
h q24[1];
cx q24[0],q24[1];
x q24[0];
h q24[0];
y q24[1];
cx q24[0],q24[1];
sdg q24[0];
h q24[0];
s q24[0];
h q24[1];
s q24[1];
z q24[1];
sdg q24[1];
h q24[1];
cx q24[0],q24[1];
x q24[0];
s q24[0];
h q24[0];
s q24[0];
h q24[1];
s q24[1];
x q24[1];
s q24[1];
h q24[1];
s q24[1];
cx q24[0],q24[1];
sdg q24[0];
h q24[0];
sdg q24[1];
h q24[1];
cx q24[0],q24[1];
h q24[0];
s q24[0];
h q24[1];
s q24[1];
cx q24[0],q24[1];
y q24[0];
sdg q24[0];
h q24[0];
x q24[1];
s q24[1];
h q24[1];
s q24[1];
cx q24[0],q24[1];
sdg q24[0];
h q24[0];
h q24[1];
s q24[1];
cx q24[0],q24[1];
h q24[0];
s q24[0];
x q24[0];
z q24[1];
s q24[1];
h q24[1];
s q24[1];
cx q24[0],q24[1];
h q24[0];
s q24[0];
sdg q24[0];
h q24[0];
sdg q24[1];
h q24[1];
z q24[1];
cx q24[0],q24[1];
y q24[0];
sdg q24[0];
h q24[0];
y q24[1];
cx q24[0],q24[1];
sdg q24[0];
h q24[0];
h q24[1];
s q24[1];
cx q24[0],q24[1];
h q24[0];
s q24[0];
cx q24[0],q24[1];
sdg q24[0];
h q24[0];
h q24[1];
s q24[1];
cx q24[0],q24[1];
h q24[0];
s q24[0];
y q24[0];
sdg q24[0];
h q24[0];
x q24[1];
cx q24[0],q24[1];
sdg q24[0];
h q24[0];
sdg q24[1];
h q24[1];
cx q24[0],q24[1];
h q24[0];
s q24[0];
h q24[0];
s q24[0];
h q24[1];
s q24[1];
z q24[1];
cx q24[0],q24[1];
h q24[0];
s q24[0];
y q24[0];
z q24[0];
h q24[0];
s q24[0];
sdg q24[1];
h q24[1];
s q24[1];
h q24[1];
s q24[1];
h q24[1];
cx q24[0],q24[1];
x q24[0];
s q24[0];
h q24[1];
s q24[1];
x q24[1];
sdg q24[1];
h q24[1];
cx q24[0],q24[1];
h q24[0];
s q24[0];
sdg q24[1];
h q24[1];
cx q24[0],q24[1];
z q24[0];
h q24[0];
h q24[1];
s q24[1];
cx q24[0],q24[1];
h q24[0];
s q24[0];
z q24[0];
s q24[0];
z q24[1];
h q24[1];
cx q24[0],q24[1];
h q24[0];
s q24[0];
sdg q24[1];
h q24[1];
cx q24[0],q24[1];
x q24[0];
s q24[0];
z q24[1];
h q24[1];
cx q24[0],q24[1];
sdg q24[0];
h q24[0];
sdg q24[1];
h q24[1];
cx q24[0],q24[1];
h q24[0];
s q24[0];
s q24[0];
h q24[1];
s q24[1];
z q24[1];
h q24[1];
cx q24[0],q24[1];
h q24[0];
s q24[0];
sdg q24[1];
h q24[1];
cx q24[0],q24[1];
z q24[0];
z q24[1];
s q24[1];
h q24[1];
s q24[1];
cx q24[0],q24[1];
h q24[0];
s q24[0];
z q24[0];
sdg q24[0];
h q24[0];
x q24[1];
sdg q24[1];
h q24[1];
cx q24[0],q24[1];
sdg q24[0];
h q24[0];
y q24[0];
h q24[0];
sdg q24[1];
h q24[1];
y q24[1];
cx q24[0],q24[1];
sdg q24[0];
h q24[0];
y q24[0];
sdg q24[0];
h q24[0];
sdg q24[1];
h q24[1];
cx q24[0],q24[1];
h q24[0];
s q24[0];
sdg q24[1];
h q24[1];
cx q24[0],q24[1];
y q24[0];
h q24[0];
h q24[1];
s q24[1];
x q24[1];
cx q24[0],q24[1];
x q24[0];
s q24[0];
sdg q24[1];
h q24[1];
z q24[1];
h q24[1];
s q24[1];
cx q24[0],q24[1];
sdg q24[0];
h q24[0];
h q24[1];
s q24[1];
cx q24[0],q24[1];
h q24[0];
s q24[0];
y q24[0];
h q24[0];
x q24[1];
s q24[1];
h q24[1];
s q24[1];
cx q24[0],q24[1];
h q24[0];
s q24[0];
h q24[1];
s q24[1];
cx q24[0],q24[1];
z q24[0];
s q24[0];
z q24[1];
s q24[1];
cx q24[0],q24[1];
h q24[0];
s q24[0];
sdg q24[1];
h q24[1];
cx q24[0],q24[1];
h q24[0];
s q24[0];
h q24[1];
s q24[1];
h q24[1];
cx q24[0],q24[1];
y q24[0];
sdg q24[1];
h q24[1];
z q24[1];
s q24[1];
h q24[1];
s q24[1];
cx q24[0],q24[1];
sdg q24[0];
h q24[0];
h q24[1];
s q24[1];
cx q24[0],q24[1];
x q24[0];
h q24[0];
s q24[0];
x q24[1];
s q24[1];
cx q24[0],q24[1];
sdg q24[0];
h q24[0];
sdg q24[1];
h q24[1];
cx q24[0],q24[1];
h q24[0];
s q24[0];
z q24[0];
s q24[0];
h q24[1];
s q24[1];
sdg q24[1];
h q24[1];
cx q24[0],q24[1];
sdg q24[0];
h q24[0];
h q24[1];
s q24[1];
cx q24[0],q24[1];
h q24[0];
s q24[0];
x q24[0];
sdg q24[0];
h q24[0];
z q24[1];
s q24[1];
h q24[1];
s q24[1];
cx q24[0],q24[1];
sdg q24[0];
h q24[0];
sdg q24[1];
h q24[1];
cx q24[0],q24[1];
y q24[0];
s q24[0];
h q24[0];
s q24[0];
h q24[1];
s q24[1];
y q24[1];
s q24[1];
cx q24[0],q24[1];
x q24[0];
s q24[0];
h q24[1];
s q24[1];
y q24[1];
cx q24[0],q24[1];
h q24[1];
s q24[1];
