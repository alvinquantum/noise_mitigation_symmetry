OPENQASM 2.0;
include "qelib1.inc";
qreg q19[2];
y q19[0];
sdg q19[0];
h q19[0];
x q19[1];
s q19[1];
cx q19[0],q19[1];
sdg q19[0];
h q19[0];
h q19[1];
s q19[1];
cx q19[0],q19[1];
h q19[0];
s q19[0];
z q19[0];
sdg q19[0];
h q19[0];
z q19[1];
h q19[1];
cx q19[0],q19[1];
sdg q19[0];
h q19[0];
sdg q19[1];
h q19[1];
cx q19[0],q19[1];
h q19[0];
s q19[0];
y q19[0];
s q19[0];
h q19[0];
s q19[0];
h q19[1];
s q19[1];
z q19[1];
s q19[1];
cx q19[0],q19[1];
sdg q19[0];
h q19[0];
sdg q19[1];
h q19[1];
cx q19[0],q19[1];
x q19[0];
h q19[0];
h q19[1];
s q19[1];
y q19[1];
h q19[1];
cx q19[0],q19[1];
x q19[0];
h q19[1];
s q19[1];
cx q19[0],q19[1];
sdg q19[0];
h q19[0];
sdg q19[1];
h q19[1];
cx q19[0],q19[1];
h q19[0];
s q19[0];
x q19[0];
h q19[0];
h q19[1];
s q19[1];
y q19[1];
sdg q19[1];
h q19[1];
cx q19[0],q19[1];
h q19[0];
s q19[0];
h q19[1];
s q19[1];
cx q19[0],q19[1];
sdg q19[0];
h q19[0];
x q19[0];
s q19[0];
h q19[0];
s q19[0];
z q19[1];
s q19[1];
h q19[1];
s q19[1];
s q19[1];
cx q19[0],q19[1];
sdg q19[0];
h q19[0];
sdg q19[1];
h q19[1];
cx q19[0],q19[1];
h q19[0];
s q19[0];
z q19[0];
sdg q19[0];
h q19[0];
z q19[1];
cx q19[0],q19[1];
sdg q19[0];
h q19[0];
z q19[0];
h q19[0];
s q19[0];
h q19[1];
s q19[1];
x q19[1];
sdg q19[1];
h q19[1];
cx q19[0],q19[1];
z q19[0];
s q19[0];
y q19[1];
cx q19[0],q19[1];
sdg q19[0];
h q19[0];
h q19[1];
s q19[1];
cx q19[0],q19[1];
h q19[0];
s q19[0];
y q19[0];
h q19[0];
s q19[0];
h q19[1];
cx q19[0],q19[1];
sdg q19[0];
h q19[0];
y q19[1];
h q19[1];
s q19[1];
cx q19[0],q19[1];
sdg q19[0];
h q19[0];
z q19[0];
s q19[0];
sdg q19[1];
h q19[1];
y q19[1];
s q19[1];
h q19[1];
s q19[1];
cx q19[0],q19[1];
sdg q19[0];
h q19[0];
sdg q19[1];
h q19[1];
cx q19[0],q19[1];
x q19[0];
h q19[0];
h q19[1];
s q19[1];
z q19[1];
s q19[1];
h q19[1];
s q19[1];
cx q19[0],q19[1];
sdg q19[0];
h q19[0];
h q19[1];
s q19[1];
cx q19[0],q19[1];
y q19[0];
h q19[0];
y q19[1];
s q19[1];
h q19[1];
s q19[1];
cx q19[0],q19[1];
sdg q19[0];
h q19[0];
h q19[1];
s q19[1];
cx q19[0],q19[1];
y q19[0];
h q19[0];
y q19[1];
s q19[1];
cx q19[0],q19[1];
x q19[0];
x q19[1];
s q19[1];
cx q19[0],q19[1];
h q19[1];
s q19[1];
y q19[1];
s q19[1];
h q19[1];
s q19[1];
cx q19[0],q19[1];
h q19[0];
s q19[0];
sdg q19[1];
h q19[1];
cx q19[0],q19[1];
s q19[0];
h q19[1];
s q19[1];
x q19[1];
h q19[1];
s q19[1];
cx q19[0],q19[1];
h q19[0];
s q19[0];
sdg q19[1];
h q19[1];
cx q19[0],q19[1];
s q19[0];
x q19[1];
s q19[1];
cx q19[0],q19[1];
sdg q19[0];
h q19[0];
sdg q19[1];
h q19[1];
cx q19[0],q19[1];
y q19[0];
h q19[0];
s q19[0];
x q19[1];
h q19[1];
s q19[1];
cx q19[0],q19[1];
h q19[1];
s q19[1];
y q19[1];
s q19[1];
h q19[1];
s q19[1];
cx q19[0],q19[1];
y q19[0];
s q19[0];
h q19[1];
s q19[1];
z q19[1];
s q19[1];
h q19[1];
s q19[1];
cx q19[0],q19[1];
sdg q19[0];
h q19[0];
sdg q19[1];
h q19[1];
cx q19[0],q19[1];
s q19[0];
h q19[1];
s q19[1];
y q19[1];
sdg q19[1];
h q19[1];
cx q19[0],q19[1];
sdg q19[0];
h q19[0];
z q19[0];
s q19[0];
h q19[0];
s q19[0];
h q19[1];
s q19[1];
y q19[1];
h q19[1];
s q19[1];
cx q19[0],q19[1];
sdg q19[0];
h q19[0];
sdg q19[1];
h q19[1];
cx q19[0],q19[1];
h q19[0];
s q19[0];
x q19[1];
h q19[1];
cx q19[0],q19[1];
sdg q19[0];
h q19[0];
x q19[0];
sdg q19[0];
h q19[0];
y q19[1];
sdg q19[1];
h q19[1];
cx q19[0],q19[1];
z q19[0];
sdg q19[0];
h q19[0];
sdg q19[1];
h q19[1];
z q19[1];
s q19[1];
cx q19[0],q19[1];
sdg q19[0];
h q19[0];
sdg q19[1];
h q19[1];
cx q19[0],q19[1];
x q19[0];
s q19[0];
h q19[0];
s q19[0];
h q19[1];
s q19[1];
cx q19[0],q19[1];
y q19[0];
sdg q19[0];
h q19[0];
h q19[1];
s q19[1];
z q19[1];
cx q19[0],q19[1];
h q19[0];
s q19[0];
y q19[0];
s q19[0];
s q19[1];
h q19[1];
s q19[1];
cx q19[0],q19[1];
h q19[0];
s q19[0];
h q19[1];
s q19[1];
cx q19[0],q19[1];
z q19[0];
x q19[1];
cx q19[0],q19[1];
h q19[0];
s q19[0];
sdg q19[1];
h q19[1];
cx q19[0],q19[1];
y q19[0];
sdg q19[0];
h q19[0];
x q19[1];
s q19[1];
h q19[1];
s q19[1];
cx q19[0],q19[1];
h q19[0];
s q19[0];
h q19[1];
s q19[1];
cx q19[0],q19[1];
y q19[0];
s q19[0];
h q19[0];
s q19[0];
z q19[1];
s q19[1];
h q19[1];
s q19[1];
cx q19[0],q19[1];
h q19[0];
s q19[0];
z q19[0];
sdg q19[0];
h q19[0];
h q19[1];
cx q19[0],q19[1];
sdg q19[0];
h q19[0];
x q19[0];
h q19[1];
s q19[1];
h q19[1];
cx q19[0],q19[1];
h q19[0];
s q19[0];
h q19[1];
s q19[1];
cx q19[0],q19[1];
y q19[0];
h q19[0];
s q19[0];
cx q19[0],q19[1];
h q19[0];
s q19[0];
x q19[0];
sdg q19[1];
h q19[1];
s q19[1];
h q19[1];
s q19[1];
cx q19[0],q19[1];
x q19[0];
x q19[0];
h q19[0];
h q19[1];
s q19[1];
y q19[1];
s q19[1];
h q19[1];
s q19[1];
y q19[1];
h q19[1];
s q19[1];
cx q19[0],q19[1];
h q19[0];
s q19[0];
z q19[0];
s q19[0];
x q19[1];
s q19[1];
cx q19[0],q19[1];
h q19[0];
s q19[0];
sdg q19[1];
h q19[1];
cx q19[0],q19[1];
z q19[0];
s q19[0];
h q19[0];
s q19[0];
h q19[1];
s q19[1];
z q19[1];
s q19[1];
h q19[1];
s q19[1];
cx q19[0],q19[1];
sdg q19[0];
h q19[0];
sdg q19[1];
h q19[1];
cx q19[0],q19[1];
y q19[0];
x q19[1];
s q19[1];
h q19[1];
s q19[1];
cx q19[0],q19[1];
sdg q19[0];
h q19[0];
sdg q19[1];
h q19[1];
cx q19[0],q19[1];
y q19[0];
h q19[0];
y q19[1];
h q19[1];
cx q19[0],q19[1];
sdg q19[0];
h q19[0];
sdg q19[1];
h q19[1];
cx q19[0],q19[1];
h q19[0];
s q19[0];
x q19[0];
h q19[0];
s q19[1];
h q19[1];
s q19[1];
cx q19[0],q19[1];
sdg q19[0];
h q19[0];
sdg q19[1];
h q19[1];
cx q19[0],q19[1];
h q19[0];
s q19[0];
h q19[1];
s q19[1];
cx q19[0],q19[1];
z q19[0];
s q19[0];
h q19[0];
s q19[0];
y q19[1];
s q19[1];
cx q19[0],q19[1];
sdg q19[0];
h q19[0];
sdg q19[1];
h q19[1];
cx q19[0],q19[1];
y q19[0];
h q19[0];
s q19[0];
z q19[1];
sdg q19[1];
h q19[1];
cx q19[0],q19[1];
h q19[0];
s q19[0];
z q19[0];
s q19[0];
h q19[0];
s q19[0];
h q19[1];
s q19[1];
z q19[1];
s q19[1];
h q19[1];
s q19[1];
cx q19[0],q19[1];
sdg q19[0];
h q19[0];
sdg q19[1];
h q19[1];
cx q19[0],q19[1];
h q19[0];
s q19[0];
y q19[1];
s q19[1];
h q19[1];
s q19[1];
cx q19[0],q19[1];
sdg q19[0];
h q19[0];
h q19[1];
s q19[1];
cx q19[0],q19[1];
z q19[0];
z q19[1];
s q19[1];
cx q19[0],q19[1];
sdg q19[0];
h q19[0];
sdg q19[1];
h q19[1];
cx q19[0],q19[1];
x q19[0];
z q19[1];
sdg q19[1];
h q19[1];
cx q19[0],q19[1];
sdg q19[0];
h q19[0];
sdg q19[0];
h q19[0];
sdg q19[1];
h q19[1];
x q19[1];
s q19[1];
cx q19[0],q19[1];
h q19[1];
s q19[1];
z q19[1];
sdg q19[1];
h q19[1];
cx q19[0],q19[1];
h q19[0];
s q19[0];
z q19[0];
sdg q19[0];
h q19[0];
sdg q19[1];
h q19[1];
h q19[1];
s q19[1];
cx q19[0],q19[1];
h q19[0];
s q19[0];
s q19[0];
h q19[0];
s q19[0];
h q19[1];
s q19[1];
s q19[1];
cx q19[0],q19[1];
sdg q19[0];
h q19[0];
sdg q19[1];
h q19[1];
cx q19[0],q19[1];
h q19[0];
s q19[0];
z q19[0];
h q19[0];
x q19[1];
h q19[1];
cx q19[0],q19[1];
h q19[0];
s q19[0];
sdg q19[1];
h q19[1];
cx q19[0],q19[1];
x q19[0];
s q19[0];
h q19[1];
s q19[1];
x q19[1];
s q19[1];
h q19[1];
s q19[1];
cx q19[0],q19[1];
h q19[0];
s q19[0];
y q19[0];
s q19[0];
h q19[0];
s q19[0];
s q19[1];
cx q19[0],q19[1];
sdg q19[0];
h q19[0];
h q19[1];
s q19[1];
cx q19[0],q19[1];
s q19[0];
h q19[0];
s q19[0];
x q19[1];
cx q19[0],q19[1];
sdg q19[0];
h q19[0];
y q19[0];
sdg q19[0];
h q19[0];
sdg q19[1];
h q19[1];
z q19[1];
h q19[1];
cx q19[0],q19[1];
x q19[0];
h q19[0];
x q19[1];
h q19[1];
s q19[1];
cx q19[0],q19[1];
y q19[0];
h q19[0];
s q19[0];
h q19[1];
s q19[1];
z q19[1];
h q19[1];
s q19[1];
cx q19[0],q19[1];
sdg q19[0];
h q19[0];
sdg q19[1];
h q19[1];
cx q19[0],q19[1];
h q19[0];
s q19[0];
y q19[1];
s q19[1];
cx q19[0],q19[1];
h q19[0];
s q19[0];
sdg q19[1];
h q19[1];
cx q19[0],q19[1];
h q19[0];
z q19[1];
cx q19[0],q19[1];
h q19[0];
s q19[0];
sdg q19[1];
h q19[1];
cx q19[0],q19[1];
z q19[0];
h q19[0];
z q19[1];
h q19[1];
s q19[1];
cx q19[0],q19[1];
sdg q19[0];
h q19[0];
sdg q19[1];
h q19[1];
cx q19[0],q19[1];
h q19[0];
s q19[0];
z q19[0];
sdg q19[0];
h q19[0];
h q19[1];
s q19[1];
h q19[1];
s q19[1];
cx q19[0],q19[1];
h q19[0];
s q19[0];
sdg q19[1];
h q19[1];
cx q19[0],q19[1];
x q19[0];
s q19[0];
h q19[0];
s q19[0];
h q19[1];
s q19[1];
h q19[1];
s q19[1];
cx q19[0],q19[1];
h q19[0];
s q19[0];
y q19[0];
s q19[0];
x q19[1];
s q19[1];
h q19[1];
s q19[1];
cx q19[0],q19[1];
sdg q19[0];
h q19[0];
sdg q19[1];
h q19[1];
cx q19[0],q19[1];
h q19[0];
s q19[0];
z q19[0];
h q19[0];
s q19[0];
h q19[1];
s q19[1];
z q19[1];
sdg q19[1];
h q19[1];
cx q19[0],q19[1];
h q19[0];
s q19[0];
sdg q19[1];
h q19[1];
cx q19[0],q19[1];
h q19[0];
y q19[1];
sdg q19[1];
h q19[1];
cx q19[0],q19[1];
z q19[0];
s q19[0];
h q19[0];
s q19[0];
h q19[1];
s q19[1];
h q19[1];
s q19[1];
cx q19[0],q19[1];
h q19[0];
s q19[0];
z q19[0];
h q19[0];
sdg q19[1];
h q19[1];
x q19[1];
cx q19[0],q19[1];
h q19[0];
s q19[0];
y q19[0];
h q19[1];
s q19[1];
y q19[1];
cx q19[0],q19[1];
sdg q19[0];
h q19[0];
sdg q19[1];
h q19[1];
cx q19[0],q19[1];
z q19[0];
h q19[0];
y q19[1];
s q19[1];
h q19[1];
s q19[1];
cx q19[0],q19[1];
sdg q19[0];
h q19[0];
y q19[0];
s q19[0];
h q19[1];
s q19[1];
x q19[1];
s q19[1];
h q19[1];
s q19[1];
cx q19[0],q19[1];
x q19[0];
h q19[0];
s q19[0];
sdg q19[1];
h q19[1];
x q19[1];
cx q19[0],q19[1];
sdg q19[0];
h q19[0];
sdg q19[1];
h q19[1];
cx q19[0],q19[1];
h q19[0];
s q19[0];
h q19[1];
s q19[1];
cx q19[0],q19[1];
h q19[1];
s q19[1];
cx q19[0],q19[1];
sdg q19[0];
h q19[0];
y q19[0];
h q19[1];
s q19[1];
z q19[1];
h q19[1];
cx q19[0],q19[1];
x q19[0];
h q19[1];
s q19[1];
x q19[1];
s q19[1];
cx q19[0],q19[1];
h q19[0];
s q19[0];
sdg q19[1];
h q19[1];
cx q19[0],q19[1];
x q19[0];
h q19[0];
s q19[0];
h q19[1];
s q19[1];
x q19[1];
h q19[1];
s q19[1];
cx q19[0],q19[1];
x q19[0];
h q19[1];
s q19[1];
x q19[1];
cx q19[0],q19[1];
sdg q19[0];
h q19[0];
sdg q19[1];
h q19[1];
cx q19[0],q19[1];
y q19[0];
y q19[1];
sdg q19[1];
h q19[1];
cx q19[0],q19[1];
h q19[0];
s q19[0];
h q19[1];
s q19[1];
cx q19[0],q19[1];
x q19[0];
h q19[0];
y q19[1];
s q19[1];
cx q19[0],q19[1];
sdg q19[1];
h q19[1];
