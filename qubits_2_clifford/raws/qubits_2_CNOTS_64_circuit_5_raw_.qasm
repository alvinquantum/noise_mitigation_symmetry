OPENQASM 2.0;
include "qelib1.inc";
qreg q4[2];
z q4[1];
h q4[1];
cx q4[0],q4[1];
z q4[0];
sdg q4[0];
h q4[0];
z q4[0];
h q4[0];
s q4[0];
y q4[1];
z q4[1];
h q4[1];
s q4[1];
cx q4[0],q4[1];
sdg q4[0];
h q4[0];
sdg q4[1];
h q4[1];
cx q4[0],q4[1];
h q4[0];
s q4[0];
h q4[0];
h q4[1];
s q4[1];
y q4[1];
s q4[1];
h q4[1];
s q4[1];
cx q4[0],q4[1];
sdg q4[0];
h q4[0];
h q4[1];
s q4[1];
cx q4[0],q4[1];
x q4[0];
s q4[0];
h q4[0];
s q4[0];
x q4[1];
cx q4[0],q4[1];
h q4[0];
s q4[0];
sdg q4[1];
h q4[1];
cx q4[0],q4[1];
z q4[0];
s q4[0];
h q4[0];
s q4[0];
h q4[1];
s q4[1];
x q4[1];
s q4[1];
cx q4[0],q4[1];
sdg q4[0];
h q4[0];
sdg q4[1];
h q4[1];
cx q4[0],q4[1];
h q4[0];
z q4[1];
cx q4[0],q4[1];
h q4[0];
s q4[0];
h q4[1];
s q4[1];
cx q4[0],q4[1];
x q4[0];
sdg q4[0];
h q4[0];
x q4[1];
sdg q4[1];
h q4[1];
cx q4[0],q4[1];
sdg q4[0];
h q4[0];
sdg q4[1];
h q4[1];
cx q4[0],q4[1];
h q4[0];
s q4[0];
z q4[0];
x q4[1];
sdg q4[1];
h q4[1];
cx q4[0],q4[1];
sdg q4[0];
h q4[0];
x q4[0];
s q4[0];
h q4[0];
s q4[0];
x q4[1];
s q4[1];
h q4[1];
s q4[1];
cx q4[0],q4[1];
sdg q4[0];
h q4[0];
x q4[0];
h q4[1];
s q4[1];
cx q4[0],q4[1];
sdg q4[0];
h q4[0];
sdg q4[1];
h q4[1];
cx q4[0],q4[1];
h q4[0];
s q4[0];
y q4[0];
sdg q4[0];
h q4[0];
h q4[1];
s q4[1];
s q4[1];
h q4[1];
s q4[1];
cx q4[0],q4[1];
sdg q4[0];
h q4[0];
sdg q4[1];
h q4[1];
cx q4[0],q4[1];
h q4[0];
h q4[1];
s q4[1];
x q4[1];
h q4[1];
s q4[1];
cx q4[0],q4[1];
h q4[0];
s q4[0];
z q4[0];
s q4[0];
y q4[0];
h q4[0];
h q4[1];
s q4[1];
z q4[1];
z q4[1];
sdg q4[1];
h q4[1];
cx q4[0],q4[1];
h q4[0];
s q4[0];
sdg q4[1];
h q4[1];
cx q4[0],q4[1];
x q4[1];
cx q4[0],q4[1];
x q4[0];
s q4[0];
h q4[1];
s q4[1];
cx q4[0],q4[1];
sdg q4[0];
h q4[0];
h q4[1];
s q4[1];
cx q4[0],q4[1];
h q4[0];
s q4[0];
z q4[0];
s q4[0];
z q4[1];
cx q4[0],q4[1];
sdg q4[0];
h q4[0];
h q4[1];
s q4[1];
cx q4[0],q4[1];
h q4[0];
s q4[0];
y q4[0];
h q4[0];
z q4[1];
s q4[1];
cx q4[0],q4[1];
h q4[0];
s q4[0];
x q4[0];
s q4[0];
y q4[0];
h q4[1];
s q4[1];
h q4[1];
cx q4[0],q4[1];
sdg q4[0];
h q4[0];
sdg q4[1];
h q4[1];
cx q4[0],q4[1];
h q4[0];
s q4[0];
h q4[1];
s q4[1];
cx q4[0],q4[1];
y q4[0];
s q4[0];
h q4[0];
s q4[0];
x q4[1];
s q4[1];
cx q4[0],q4[1];
y q4[0];
h q4[0];
s q4[0];
x q4[0];
s q4[0];
s q4[1];
h q4[1];
s q4[1];
y q4[1];
h q4[1];
cx q4[0],q4[1];
sdg q4[0];
h q4[0];
sdg q4[1];
h q4[1];
cx q4[0],q4[1];
h q4[0];
s q4[0];
h q4[1];
s q4[1];
cx q4[0],q4[1];
s q4[0];
z q4[1];
h q4[1];
cx q4[0],q4[1];
sdg q4[0];
h q4[0];
sdg q4[1];
h q4[1];
cx q4[0],q4[1];
h q4[0];
s q4[0];
h q4[0];
s q4[0];
h q4[1];
s q4[1];
x q4[1];
cx q4[0],q4[1];
sdg q4[0];
h q4[0];
y q4[0];
h q4[0];
z q4[1];
s q4[1];
cx q4[0],q4[1];
y q4[0];
h q4[0];
s q4[0];
h q4[1];
s q4[1];
cx q4[0],q4[1];
sdg q4[0];
h q4[0];
sdg q4[1];
h q4[1];
cx q4[0],q4[1];
h q4[0];
s q4[0];
y q4[0];
h q4[1];
s q4[1];
x q4[1];
s q4[1];
cx q4[0],q4[1];
sdg q4[0];
h q4[0];
z q4[0];
h q4[0];
h q4[1];
s q4[1];
y q4[1];
h q4[1];
s q4[1];
cx q4[0],q4[1];
sdg q4[0];
h q4[0];
sdg q4[1];
h q4[1];
cx q4[0],q4[1];
h q4[0];
s q4[0];
z q4[0];
h q4[1];
s q4[1];
z q4[1];
cx q4[0],q4[1];
h q4[0];
s q4[0];
s q4[0];
h q4[0];
s q4[0];
sdg q4[1];
h q4[1];
h q4[1];
s q4[1];
cx q4[0],q4[1];
sdg q4[0];
h q4[0];
x q4[0];
h q4[0];
h q4[1];
s q4[1];
x q4[1];
sdg q4[1];
h q4[1];
cx q4[0],q4[1];
sdg q4[0];
h q4[0];
sdg q4[1];
h q4[1];
cx q4[0],q4[1];
h q4[0];
s q4[0];
x q4[0];
s q4[0];
h q4[0];
s q4[0];
h q4[1];
cx q4[0],q4[1];
sdg q4[0];
h q4[0];
sdg q4[1];
h q4[1];
cx q4[0],q4[1];
h q4[0];
s q4[0];
x q4[0];
s q4[0];
h q4[1];
s q4[1];
x q4[1];
s q4[1];
cx q4[0],q4[1];
x q4[0];
s q4[0];
h q4[0];
s q4[0];
sdg q4[1];
h q4[1];
z q4[1];
s q4[1];
h q4[1];
s q4[1];
cx q4[0],q4[1];
h q4[0];
s q4[0];
sdg q4[1];
h q4[1];
s q4[1];
cx q4[0],q4[1];
h q4[0];
s q4[0];
z q4[0];
s q4[0];
h q4[0];
s q4[0];
x q4[1];
s q4[1];
h q4[1];
s q4[1];
cx q4[0],q4[1];
sdg q4[0];
h q4[0];
sdg q4[1];
h q4[1];
cx q4[0],q4[1];
h q4[0];
z q4[1];
cx q4[0],q4[1];
sdg q4[0];
h q4[0];
sdg q4[1];
h q4[1];
cx q4[0],q4[1];
h q4[0];
h q4[1];
s q4[1];
x q4[1];
h q4[1];
cx q4[0],q4[1];
h q4[0];
s q4[0];
sdg q4[1];
h q4[1];
cx q4[0],q4[1];
x q4[0];
h q4[1];
s q4[1];
y q4[1];
sdg q4[1];
h q4[1];
cx q4[0],q4[1];
sdg q4[0];
h q4[0];
h q4[1];
s q4[1];
cx q4[0],q4[1];
h q4[0];
s q4[0];
z q4[1];
cx q4[0],q4[1];
h q4[0];
s q4[0];
sdg q4[1];
h q4[1];
cx q4[0],q4[1];
h q4[0];
s q4[0];
h q4[1];
s q4[1];
y q4[1];
sdg q4[1];
h q4[1];
cx q4[0],q4[1];
sdg q4[0];
h q4[0];
sdg q4[1];
h q4[1];
