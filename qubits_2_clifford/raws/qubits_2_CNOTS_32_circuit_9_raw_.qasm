OPENQASM 2.0;
include "qelib1.inc";
qreg q8[2];
x q8[0];
s q8[0];
h q8[0];
s q8[0];
z q8[1];
h q8[1];
s q8[1];
cx q8[0],q8[1];
h q8[0];
s q8[0];
y q8[0];
sdg q8[1];
h q8[1];
sdg q8[1];
h q8[1];
cx q8[0],q8[1];
sdg q8[0];
h q8[0];
y q8[0];
s q8[0];
h q8[1];
s q8[1];
cx q8[0],q8[1];
sdg q8[0];
h q8[0];
sdg q8[1];
h q8[1];
cx q8[0],q8[1];
h q8[0];
s q8[0];
h q8[1];
s q8[1];
cx q8[0],q8[1];
z q8[0];
sdg q8[0];
h q8[0];
x q8[1];
h q8[1];
cx q8[0],q8[1];
sdg q8[0];
h q8[0];
z q8[0];
sdg q8[0];
h q8[0];
x q8[1];
s q8[1];
h q8[1];
s q8[1];
cx q8[0],q8[1];
sdg q8[1];
h q8[1];
y q8[1];
h q8[1];
s q8[1];
cx q8[0],q8[1];
h q8[0];
s q8[0];
z q8[0];
sdg q8[1];
h q8[1];
z q8[1];
sdg q8[1];
h q8[1];
cx q8[0],q8[1];
h q8[0];
s q8[0];
x q8[0];
h q8[0];
sdg q8[1];
h q8[1];
sdg q8[1];
h q8[1];
cx q8[0],q8[1];
sdg q8[0];
h q8[0];
sdg q8[1];
h q8[1];
cx q8[0],q8[1];
h q8[0];
s q8[0];
x q8[0];
h q8[0];
s q8[0];
s q8[1];
cx q8[0],q8[1];
h q8[0];
s q8[0];
y q8[0];
h q8[1];
s q8[1];
x q8[1];
sdg q8[1];
h q8[1];
cx q8[0],q8[1];
sdg q8[0];
h q8[0];
s q8[0];
h q8[0];
s q8[0];
cx q8[0],q8[1];
y q8[0];
h q8[0];
z q8[1];
s q8[1];
cx q8[0],q8[1];
h q8[0];
s q8[0];
y q8[0];
s q8[0];
z q8[1];
h q8[1];
s q8[1];
cx q8[0],q8[1];
z q8[0];
s q8[0];
h q8[0];
s q8[0];
y q8[1];
s q8[1];
h q8[1];
s q8[1];
cx q8[0],q8[1];
sdg q8[0];
h q8[0];
sdg q8[1];
h q8[1];
cx q8[0],q8[1];
h q8[0];
s q8[0];
z q8[0];
sdg q8[0];
h q8[0];
h q8[1];
s q8[1];
s q8[1];
cx q8[0],q8[1];
sdg q8[0];
h q8[0];
sdg q8[1];
h q8[1];
cx q8[0],q8[1];
x q8[0];
s q8[0];
h q8[1];
s q8[1];
x q8[1];
sdg q8[1];
h q8[1];
cx q8[0],q8[1];
x q8[0];
s q8[0];
h q8[0];
s q8[0];
sdg q8[1];
h q8[1];
z q8[1];
h q8[1];
s q8[1];
cx q8[0],q8[1];
sdg q8[0];
h q8[0];
sdg q8[1];
h q8[1];
cx q8[0],q8[1];
h q8[0];
s q8[0];
y q8[0];
sdg q8[0];
h q8[0];
y q8[1];
s q8[1];
h q8[1];
s q8[1];
cx q8[0],q8[1];
sdg q8[0];
h q8[0];
h q8[1];
s q8[1];
cx q8[0],q8[1];
y q8[0];
s q8[0];
sdg q8[1];
h q8[1];
cx q8[0],q8[1];
sdg q8[0];
h q8[0];
x q8[0];
h q8[0];
s q8[0];
h q8[1];
s q8[1];
x q8[1];
s q8[1];
cx q8[0],q8[1];
y q8[1];
h q8[1];
cx q8[0],q8[1];
h q8[0];
s q8[0];
y q8[0];
s q8[0];
h q8[0];
s q8[0];
y q8[1];
h q8[1];
cx q8[0],q8[1];
y q8[0];
s q8[0];
z q8[1];
s q8[1];
cx q8[0],q8[1];
sdg q8[0];
h q8[0];
h q8[1];
s q8[1];
cx q8[0],q8[1];
z q8[0];
sdg q8[0];
h q8[0];
s q8[0];
y q8[1];
h q8[1];
x q8[1];
h q8[1];
s q8[1];
cx q8[0],q8[1];
h q8[1];
s q8[1];
