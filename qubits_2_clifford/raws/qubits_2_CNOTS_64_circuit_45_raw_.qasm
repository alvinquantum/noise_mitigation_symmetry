OPENQASM 2.0;
include "qelib1.inc";
qreg q44[2];
z q44[1];
cx q44[0],q44[1];
h q44[0];
sdg q44[1];
h q44[1];
sdg q44[1];
h q44[1];
cx q44[0],q44[1];
h q44[0];
s q44[0];
z q44[0];
h q44[0];
s q44[0];
h q44[1];
s q44[1];
y q44[1];
cx q44[0],q44[1];
sdg q44[0];
h q44[0];
h q44[1];
s q44[1];
cx q44[0],q44[1];
h q44[0];
s q44[0];
z q44[0];
s q44[0];
cx q44[0],q44[1];
sdg q44[0];
h q44[0];
y q44[0];
y q44[1];
sdg q44[1];
h q44[1];
cx q44[0],q44[1];
x q44[0];
s q44[0];
h q44[0];
s q44[0];
h q44[1];
s q44[1];
x q44[1];
h q44[1];
s q44[1];
cx q44[0],q44[1];
h q44[0];
s q44[0];
h q44[1];
s q44[1];
cx q44[0],q44[1];
y q44[0];
s q44[0];
z q44[1];
s q44[1];
cx q44[0],q44[1];
sdg q44[0];
h q44[0];
sdg q44[1];
h q44[1];
cx q44[0],q44[1];
h q44[0];
s q44[0];
z q44[0];
h q44[1];
cx q44[0],q44[1];
h q44[0];
s q44[0];
sdg q44[1];
h q44[1];
cx q44[0],q44[1];
z q44[0];
s q44[0];
h q44[1];
cx q44[0],q44[1];
sdg q44[0];
h q44[0];
z q44[0];
h q44[0];
h q44[1];
s q44[1];
z q44[1];
cx q44[0],q44[1];
sdg q44[0];
h q44[0];
y q44[0];
s q44[0];
h q44[1];
s q44[1];
y q44[1];
s q44[1];
h q44[1];
s q44[1];
cx q44[0],q44[1];
sdg q44[0];
h q44[0];
h q44[1];
s q44[1];
cx q44[0],q44[1];
x q44[0];
h q44[0];
s q44[0];
x q44[1];
h q44[1];
s q44[1];
cx q44[0],q44[1];
h q44[0];
s q44[0];
y q44[0];
s q44[0];
h q44[0];
s q44[0];
x q44[0];
s q44[0];
h q44[1];
s q44[1];
h q44[1];
s q44[1];
h q44[1];
s q44[1];
cx q44[0],q44[1];
h q44[0];
s q44[0];
h q44[1];
s q44[1];
cx q44[0],q44[1];
x q44[0];
s q44[0];
x q44[1];
h q44[1];
s q44[1];
cx q44[0],q44[1];
h q44[0];
s q44[0];
z q44[0];
h q44[0];
z q44[1];
s q44[1];
cx q44[0],q44[1];
sdg q44[0];
h q44[0];
x q44[0];
s q44[0];
h q44[1];
s q44[1];
z q44[1];
sdg q44[1];
h q44[1];
cx q44[0],q44[1];
sdg q44[0];
h q44[0];
sdg q44[1];
h q44[1];
cx q44[0],q44[1];
h q44[0];
s q44[0];
x q44[0];
sdg q44[0];
h q44[0];
z q44[1];
sdg q44[1];
h q44[1];
cx q44[0],q44[1];
h q44[0];
s q44[0];
z q44[0];
s q44[0];
h q44[0];
s q44[0];
h q44[1];
s q44[1];
x q44[1];
h q44[1];
cx q44[0],q44[1];
sdg q44[0];
h q44[0];
sdg q44[1];
h q44[1];
cx q44[0],q44[1];
y q44[0];
sdg q44[0];
h q44[0];
h q44[1];
s q44[1];
x q44[1];
h q44[1];
cx q44[0],q44[1];
sdg q44[0];
h q44[0];
sdg q44[1];
h q44[1];
cx q44[0],q44[1];
x q44[0];
h q44[0];
z q44[1];
s q44[1];
cx q44[0],q44[1];
s q44[0];
h q44[0];
s q44[0];
sdg q44[1];
h q44[1];
y q44[1];
s q44[1];
h q44[1];
s q44[1];
cx q44[0],q44[1];
sdg q44[0];
h q44[0];
sdg q44[1];
h q44[1];
cx q44[0],q44[1];
h q44[0];
s q44[0];
h q44[1];
s q44[1];
cx q44[0],q44[1];
y q44[0];
s q44[0];
h q44[1];
s q44[1];
cx q44[0],q44[1];
h q44[0];
s q44[0];
h q44[1];
s q44[1];
cx q44[0],q44[1];
s q44[0];
h q44[0];
s q44[0];
y q44[1];
s q44[1];
cx q44[0],q44[1];
sdg q44[0];
h q44[0];
sdg q44[1];
h q44[1];
cx q44[0],q44[1];
h q44[0];
s q44[0];
h q44[0];
s q44[0];
h q44[1];
s q44[1];
y q44[1];
h q44[1];
s q44[1];
cx q44[0],q44[1];
sdg q44[0];
h q44[0];
sdg q44[1];
h q44[1];
cx q44[0],q44[1];
y q44[0];
s q44[0];
h q44[0];
s q44[0];
y q44[1];
h q44[1];
s q44[1];
cx q44[0],q44[1];
h q44[0];
s q44[0];
h q44[1];
s q44[1];
cx q44[0],q44[1];
z q44[0];
y q44[1];
h q44[1];
cx q44[0],q44[1];
sdg q44[0];
h q44[0];
sdg q44[1];
h q44[1];
cx q44[0],q44[1];
z q44[0];
s q44[0];
h q44[0];
s q44[0];
h q44[1];
s q44[1];
y q44[1];
s q44[1];
h q44[1];
s q44[1];
cx q44[0],q44[1];
sdg q44[0];
h q44[0];
sdg q44[1];
h q44[1];
cx q44[0],q44[1];
y q44[0];
y q44[1];
s q44[1];
cx q44[0],q44[1];
y q44[0];
s q44[0];
h q44[0];
s q44[0];
z q44[1];
s q44[1];
cx q44[0],q44[1];
sdg q44[0];
h q44[0];
sdg q44[1];
h q44[1];
cx q44[0],q44[1];
h q44[0];
y q44[1];
cx q44[0],q44[1];
h q44[0];
s q44[0];
y q44[0];
sdg q44[0];
h q44[0];
h q44[1];
s q44[1];
x q44[1];
s q44[1];
h q44[1];
s q44[1];
cx q44[0],q44[1];
sdg q44[0];
h q44[0];
z q44[0];
s q44[0];
sdg q44[1];
h q44[1];
z q44[1];
h q44[1];
cx q44[0],q44[1];
h q44[0];
s q44[0];
sdg q44[1];
h q44[1];
cx q44[0],q44[1];
s q44[0];
z q44[1];
h q44[1];
s q44[1];
cx q44[0],q44[1];
sdg q44[0];
h q44[0];
sdg q44[1];
h q44[1];
cx q44[0],q44[1];
h q44[0];
s q44[0];
x q44[0];
y q44[1];
cx q44[0],q44[1];
sdg q44[0];
h q44[0];
h q44[1];
s q44[1];
cx q44[0],q44[1];
y q44[0];
sdg q44[0];
h q44[0];
y q44[1];
s q44[1];
h q44[1];
s q44[1];
cx q44[0],q44[1];
sdg q44[0];
h q44[0];
s q44[0];
sdg q44[1];
h q44[1];
x q44[1];
h q44[1];
s q44[1];
cx q44[0],q44[1];
h q44[0];
s q44[0];
y q44[0];
h q44[1];
s q44[1];
z q44[1];
sdg q44[1];
h q44[1];
cx q44[0],q44[1];
sdg q44[0];
h q44[0];
sdg q44[1];
h q44[1];
cx q44[0],q44[1];
h q44[0];
s q44[0];
y q44[0];
s q44[0];
h q44[1];
s q44[1];
y q44[1];
h q44[1];
cx q44[0],q44[1];
h q44[0];
s q44[0];
x q44[0];
sdg q44[0];
h q44[0];
sdg q44[0];
h q44[0];
z q44[1];
s q44[1];
h q44[1];
s q44[1];
x q44[1];
s q44[1];
cx q44[0],q44[1];
h q44[0];
s q44[0];
h q44[1];
s q44[1];
cx q44[0],q44[1];
z q44[0];
h q44[0];
s q44[1];
h q44[1];
s q44[1];
cx q44[0],q44[1];
z q44[0];
sdg q44[0];
h q44[0];
y q44[1];
s q44[1];
h q44[1];
s q44[1];
cx q44[0],q44[1];
h q44[0];
s q44[0];
sdg q44[1];
h q44[1];
