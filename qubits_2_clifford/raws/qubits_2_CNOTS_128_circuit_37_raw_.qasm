OPENQASM 2.0;
include "qelib1.inc";
qreg q36[2];
x q36[0];
h q36[0];
s q36[0];
z q36[1];
h q36[1];
s q36[1];
cx q36[0],q36[1];
h q36[0];
s q36[0];
sdg q36[1];
h q36[1];
cx q36[0],q36[1];
y q36[0];
h q36[1];
s q36[1];
z q36[1];
h q36[1];
cx q36[0],q36[1];
h q36[0];
s q36[0];
z q36[0];
s q36[0];
h q36[0];
s q36[0];
sdg q36[1];
h q36[1];
h q36[1];
s q36[1];
cx q36[0],q36[1];
sdg q36[0];
h q36[0];
sdg q36[1];
h q36[1];
cx q36[0],q36[1];
h q36[0];
s q36[0];
x q36[0];
sdg q36[0];
h q36[0];
y q36[0];
h q36[0];
y q36[1];
sdg q36[1];
h q36[1];
y q36[1];
s q36[1];
h q36[1];
s q36[1];
cx q36[0],q36[1];
sdg q36[0];
h q36[0];
sdg q36[1];
h q36[1];
cx q36[0],q36[1];
h q36[0];
s q36[0];
h q36[1];
s q36[1];
cx q36[0],q36[1];
h q36[0];
z q36[1];
h q36[1];
cx q36[0],q36[1];
sdg q36[0];
h q36[0];
sdg q36[1];
h q36[1];
cx q36[0],q36[1];
x q36[0];
sdg q36[0];
h q36[0];
h q36[1];
s q36[1];
y q36[1];
s q36[1];
cx q36[0],q36[1];
sdg q36[0];
h q36[0];
h q36[1];
s q36[1];
cx q36[0],q36[1];
h q36[0];
s q36[0];
y q36[0];
x q36[1];
h q36[1];
s q36[1];
cx q36[0],q36[1];
sdg q36[0];
h q36[0];
z q36[0];
z q36[1];
h q36[1];
cx q36[0],q36[1];
sdg q36[0];
h q36[0];
y q36[0];
s q36[0];
h q36[1];
s q36[1];
y q36[1];
sdg q36[1];
h q36[1];
cx q36[0],q36[1];
h q36[0];
s q36[0];
z q36[0];
z q36[1];
sdg q36[1];
h q36[1];
cx q36[0],q36[1];
y q36[0];
h q36[0];
s q36[0];
y q36[1];
h q36[1];
s q36[1];
cx q36[0],q36[1];
sdg q36[0];
h q36[0];
h q36[1];
s q36[1];
z q36[1];
h q36[1];
cx q36[0],q36[1];
h q36[0];
s q36[0];
y q36[0];
s q36[0];
h q36[1];
s q36[1];
sdg q36[1];
h q36[1];
cx q36[0],q36[1];
h q36[0];
s q36[0];
z q36[0];
h q36[0];
sdg q36[1];
h q36[1];
x q36[1];
s q36[1];
h q36[1];
s q36[1];
cx q36[0],q36[1];
x q36[0];
cx q36[0],q36[1];
sdg q36[0];
h q36[0];
sdg q36[1];
h q36[1];
cx q36[0],q36[1];
sdg q36[0];
h q36[0];
h q36[1];
s q36[1];
x q36[1];
h q36[1];
s q36[1];
cx q36[0],q36[1];
sdg q36[0];
h q36[0];
sdg q36[1];
h q36[1];
cx q36[0],q36[1];
x q36[0];
h q36[0];
x q36[1];
s q36[1];
cx q36[0],q36[1];
sdg q36[0];
h q36[0];
s q36[0];
h q36[1];
s q36[1];
z q36[1];
sdg q36[1];
h q36[1];
cx q36[0],q36[1];
h q36[0];
s q36[0];
h q36[1];
s q36[1];
cx q36[0],q36[1];
y q36[0];
s q36[0];
h q36[0];
s q36[0];
z q36[1];
s q36[1];
cx q36[0],q36[1];
h q36[0];
s q36[0];
x q36[0];
h q36[0];
sdg q36[1];
h q36[1];
x q36[1];
sdg q36[1];
h q36[1];
cx q36[0],q36[1];
h q36[0];
s q36[0];
h q36[1];
s q36[1];
cx q36[0],q36[1];
x q36[0];
sdg q36[0];
h q36[0];
z q36[1];
cx q36[0],q36[1];
sdg q36[0];
h q36[0];
h q36[1];
s q36[1];
h q36[1];
cx q36[0],q36[1];
h q36[0];
s q36[0];
h q36[1];
s q36[1];
cx q36[0],q36[1];
z q36[0];
sdg q36[0];
h q36[0];
z q36[1];
cx q36[0],q36[1];
h q36[0];
s q36[0];
sdg q36[1];
h q36[1];
cx q36[0],q36[1];
sdg q36[0];
h q36[0];
y q36[1];
s q36[1];
cx q36[0],q36[1];
h q36[0];
s q36[0];
sdg q36[1];
h q36[1];
cx q36[0],q36[1];
x q36[0];
sdg q36[0];
h q36[0];
x q36[1];
h q36[1];
s q36[1];
cx q36[0],q36[1];
s q36[0];
h q36[1];
s q36[1];
y q36[1];
sdg q36[1];
h q36[1];
cx q36[0],q36[1];
sdg q36[0];
h q36[0];
sdg q36[1];
h q36[1];
cx q36[0],q36[1];
h q36[0];
s q36[0];
z q36[0];
h q36[1];
s q36[1];
y q36[1];
h q36[1];
cx q36[0],q36[1];
h q36[0];
s q36[0];
y q36[0];
s q36[0];
h q36[0];
s q36[0];
h q36[1];
s q36[1];
y q36[1];
cx q36[0],q36[1];
sdg q36[0];
h q36[0];
sdg q36[1];
h q36[1];
cx q36[0],q36[1];
h q36[0];
s q36[0];
x q36[0];
y q36[1];
s q36[1];
h q36[1];
s q36[1];
cx q36[0],q36[1];
h q36[0];
s q36[0];
z q36[0];
sdg q36[1];
h q36[1];
y q36[1];
sdg q36[1];
h q36[1];
cx q36[0],q36[1];
sdg q36[0];
h q36[0];
sdg q36[1];
h q36[1];
cx q36[0],q36[1];
h q36[0];
s q36[0];
h q36[1];
s q36[1];
cx q36[0],q36[1];
y q36[0];
sdg q36[0];
h q36[0];
x q36[1];
s q36[1];
h q36[1];
s q36[1];
cx q36[0],q36[1];
sdg q36[0];
h q36[0];
s q36[0];
h q36[0];
s q36[0];
sdg q36[1];
h q36[1];
z q36[1];
cx q36[0],q36[1];
h q36[0];
s q36[0];
sdg q36[1];
h q36[1];
cx q36[0],q36[1];
y q36[0];
sdg q36[0];
h q36[0];
z q36[1];
s q36[1];
h q36[1];
s q36[1];
cx q36[0],q36[1];
h q36[0];
s q36[0];
sdg q36[1];
h q36[1];
cx q36[0],q36[1];
z q36[0];
h q36[0];
h q36[1];
s q36[1];
z q36[1];
s q36[1];
cx q36[0],q36[1];
sdg q36[0];
h q36[0];
sdg q36[1];
h q36[1];
cx q36[0],q36[1];
h q36[0];
s q36[0];
y q36[0];
sdg q36[0];
h q36[0];
h q36[1];
s q36[1];
y q36[1];
h q36[1];
s q36[1];
cx q36[0],q36[1];
sdg q36[0];
h q36[0];
sdg q36[1];
h q36[1];
cx q36[0],q36[1];
x q36[0];
s q36[0];
sdg q36[0];
h q36[0];
z q36[1];
h q36[1];
y q36[1];
sdg q36[1];
h q36[1];
cx q36[0],q36[1];
z q36[0];
sdg q36[0];
h q36[0];
z q36[1];
h q36[1];
s q36[1];
cx q36[0],q36[1];
sdg q36[0];
h q36[0];
sdg q36[1];
h q36[1];
cx q36[0],q36[1];
h q36[0];
s q36[0];
x q36[0];
s q36[0];
y q36[1];
s q36[1];
h q36[1];
s q36[1];
cx q36[0],q36[1];
h q36[0];
s q36[0];
h q36[1];
s q36[1];
cx q36[0],q36[1];
z q36[1];
s q36[1];
cx q36[0],q36[1];
sdg q36[0];
h q36[0];
sdg q36[1];
h q36[1];
cx q36[0],q36[1];
sdg q36[0];
h q36[0];
x q36[1];
s q36[1];
h q36[1];
s q36[1];
cx q36[0],q36[1];
sdg q36[0];
h q36[0];
sdg q36[1];
h q36[1];
cx q36[0],q36[1];
x q36[0];
s q36[0];
h q36[1];
s q36[1];
z q36[1];
sdg q36[1];
h q36[1];
cx q36[0],q36[1];
h q36[0];
s q36[0];
h q36[1];
s q36[1];
cx q36[0],q36[1];
sdg q36[0];
h q36[0];
sdg q36[1];
h q36[1];
cx q36[0],q36[1];
h q36[0];
s q36[0];
h q36[1];
s q36[1];
cx q36[0],q36[1];
y q36[0];
s q36[0];
x q36[0];
h q36[1];
y q36[1];
sdg q36[1];
h q36[1];
cx q36[0],q36[1];
h q36[0];
s q36[0];
z q36[0];
s q36[0];
z q36[1];
s q36[1];
cx q36[0],q36[1];
sdg q36[0];
h q36[0];
h q36[1];
s q36[1];
cx q36[0],q36[1];
y q36[0];
s q36[0];
x q36[1];
s q36[1];
h q36[1];
s q36[1];
cx q36[0],q36[1];
y q36[0];
sdg q36[0];
h q36[0];
sdg q36[1];
h q36[1];
x q36[1];
s q36[1];
cx q36[0],q36[1];
h q36[0];
s q36[0];
sdg q36[0];
h q36[0];
sdg q36[1];
h q36[1];
x q36[1];
s q36[1];
cx q36[0],q36[1];
sdg q36[0];
h q36[0];
h q36[1];
s q36[1];
cx q36[0],q36[1];
x q36[0];
h q36[0];
z q36[0];
y q36[1];
sdg q36[1];
h q36[1];
z q36[1];
h q36[1];
s q36[1];
cx q36[0],q36[1];
h q36[0];
s q36[0];
sdg q36[1];
h q36[1];
cx q36[0],q36[1];
y q36[0];
h q36[0];
z q36[1];
h q36[1];
s q36[1];
cx q36[0],q36[1];
sdg q36[0];
h q36[0];
z q36[0];
s q36[0];
y q36[1];
cx q36[0],q36[1];
sdg q36[0];
h q36[0];
y q36[0];
h q36[1];
s q36[1];
z q36[1];
s q36[1];
h q36[1];
s q36[1];
cx q36[0],q36[1];
y q36[0];
s q36[0];
sdg q36[1];
h q36[1];
s q36[1];
cx q36[0],q36[1];
sdg q36[0];
h q36[0];
sdg q36[1];
h q36[1];
cx q36[0],q36[1];
y q36[0];
h q36[1];
s q36[1];
s q36[1];
h q36[1];
s q36[1];
cx q36[0],q36[1];
h q36[0];
s q36[0];
h q36[1];
s q36[1];
z q36[1];
h q36[1];
cx q36[0],q36[1];
h q36[0];
s q36[0];
sdg q36[1];
h q36[1];
cx q36[0],q36[1];
x q36[0];
h q36[1];
s q36[1];
h q36[1];
s q36[1];
cx q36[0],q36[1];
h q36[0];
s q36[0];
sdg q36[1];
h q36[1];
cx q36[0],q36[1];
z q36[0];
s q36[0];
h q36[0];
s q36[0];
h q36[1];
s q36[1];
z q36[1];
s q36[1];
h q36[1];
s q36[1];
cx q36[0],q36[1];
sdg q36[0];
h q36[0];
sdg q36[1];
h q36[1];
cx q36[0],q36[1];
h q36[0];
s q36[0];
y q36[1];
s q36[1];
cx q36[0],q36[1];
sdg q36[0];
h q36[0];
sdg q36[1];
h q36[1];
cx q36[0],q36[1];
h q36[1];
s q36[1];
y q36[1];
h q36[1];
cx q36[0],q36[1];
sdg q36[0];
h q36[0];
y q36[0];
h q36[0];
h q36[1];
s q36[1];
h q36[1];
s q36[1];
cx q36[0],q36[1];
sdg q36[0];
h q36[0];
h q36[1];
s q36[1];
cx q36[0],q36[1];
z q36[0];
s q36[1];
cx q36[0],q36[1];
sdg q36[0];
h q36[0];
s q36[0];
sdg q36[1];
h q36[1];
z q36[1];
s q36[1];
cx q36[0],q36[1];
sdg q36[0];
h q36[0];
y q36[0];
s q36[0];
s q36[1];
cx q36[0],q36[1];
z q36[0];
h q36[0];
sdg q36[1];
h q36[1];
z q36[1];
sdg q36[1];
h q36[1];
cx q36[0],q36[1];
z q36[0];
s q36[0];
h q36[0];
s q36[0];
s q36[1];
cx q36[0],q36[1];
h q36[0];
s q36[0];
sdg q36[1];
h q36[1];
cx q36[0],q36[1];
y q36[0];
h q36[0];
s q36[0];
h q36[1];
s q36[1];
x q36[1];
sdg q36[1];
h q36[1];
cx q36[0],q36[1];
sdg q36[0];
h q36[0];
h q36[1];
s q36[1];
cx q36[0],q36[1];
z q36[0];
cx q36[0],q36[1];
sdg q36[0];
h q36[0];
h q36[1];
s q36[1];
cx q36[0],q36[1];
z q36[0];
h q36[0];
sdg q36[1];
h q36[1];
cx q36[0],q36[1];
h q36[0];
s q36[0];
sdg q36[1];
h q36[1];
cx q36[0],q36[1];
z q36[0];
sdg q36[0];
h q36[0];
z q36[1];
s q36[1];
h q36[1];
s q36[1];
cx q36[0],q36[1];
h q36[0];
s q36[0];
z q36[0];
sdg q36[0];
h q36[0];
h q36[1];
s q36[1];
x q36[1];
s q36[1];
cx q36[0],q36[1];
h q36[0];
s q36[0];
h q36[0];
s q36[0];
h q36[1];
s q36[1];
y q36[1];
s q36[1];
cx q36[0],q36[1];
x q36[0];
h q36[0];
h q36[1];
s q36[1];
x q36[1];
sdg q36[1];
h q36[1];
cx q36[0],q36[1];
sdg q36[0];
h q36[0];
s q36[0];
h q36[0];
s q36[0];
y q36[1];
h q36[1];
s q36[1];
cx q36[0],q36[1];
sdg q36[0];
h q36[0];
sdg q36[1];
h q36[1];
cx q36[0],q36[1];
h q36[0];
s q36[0];
h q36[1];
s q36[1];
cx q36[0],q36[1];
s q36[0];
z q36[1];
h q36[1];
s q36[1];
cx q36[0],q36[1];
sdg q36[0];
h q36[0];
sdg q36[1];
h q36[1];
cx q36[0],q36[1];
h q36[0];
s q36[0];
z q36[0];
sdg q36[0];
h q36[0];
y q36[1];
h q36[1];
cx q36[0],q36[1];
h q36[0];
s q36[0];
z q36[0];
sdg q36[0];
h q36[0];
z q36[1];
h q36[1];
s q36[1];
cx q36[0],q36[1];
sdg q36[0];
h q36[0];
h q36[1];
s q36[1];
cx q36[0],q36[1];
z q36[0];
s q36[0];
h q36[0];
s q36[0];
z q36[1];
cx q36[0],q36[1];
h q36[0];
s q36[0];
y q36[0];
sdg q36[0];
h q36[0];
h q36[1];
s q36[1];
sdg q36[1];
h q36[1];
cx q36[0],q36[1];
h q36[1];
s q36[1];
x q36[1];
s q36[1];
h q36[1];
s q36[1];
cx q36[0],q36[1];
sdg q36[0];
h q36[0];
sdg q36[1];
h q36[1];
cx q36[0],q36[1];
s q36[0];
h q36[1];
s q36[1];
z q36[1];
h q36[1];
cx q36[0],q36[1];
s q36[0];
h q36[0];
s q36[0];
x q36[1];
s q36[1];
cx q36[0],q36[1];
sdg q36[0];
h q36[0];
y q36[0];
s q36[0];
h q36[1];
s q36[1];
s q36[1];
cx q36[0],q36[1];
sdg q36[0];
h q36[0];
sdg q36[0];
h q36[0];
h q36[1];
s q36[1];
x q36[1];
sdg q36[1];
h q36[1];
cx q36[0],q36[1];
sdg q36[0];
h q36[0];
h q36[1];
s q36[1];
cx q36[0],q36[1];
h q36[0];
s q36[0];
