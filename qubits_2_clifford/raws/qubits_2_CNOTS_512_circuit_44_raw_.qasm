OPENQASM 2.0;
include "qelib1.inc";
qreg q43[2];
z q43[0];
x q43[1];
s q43[1];
h q43[1];
s q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
h q43[1];
s q43[1];
cx q43[0],q43[1];
z q43[0];
sdg q43[0];
h q43[0];
y q43[1];
h q43[1];
s q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
h q43[1];
s q43[1];
cx q43[0],q43[1];
y q43[0];
s q43[0];
x q43[1];
s q43[1];
h q43[1];
s q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
h q43[1];
s q43[1];
cx q43[0],q43[1];
z q43[0];
h q43[0];
s q43[0];
x q43[0];
h q43[0];
x q43[1];
h q43[1];
s q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
s q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
h q43[1];
s q43[1];
cx q43[0],q43[1];
y q43[0];
h q43[0];
s q43[0];
y q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
x q43[0];
h q43[0];
z q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
y q43[0];
h q43[0];
s q43[0];
h q43[1];
s q43[1];
y q43[1];
s q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
x q43[0];
sdg q43[0];
h q43[0];
h q43[1];
s q43[1];
y q43[1];
s q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
x q43[0];
sdg q43[0];
h q43[0];
h q43[1];
s q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
z q43[0];
h q43[0];
sdg q43[1];
h q43[1];
z q43[1];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
h q43[0];
h q43[1];
s q43[1];
z q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
h q43[1];
s q43[1];
z q43[1];
s q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
z q43[0];
s q43[0];
h q43[1];
s q43[1];
x q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
s q43[0];
h q43[0];
s q43[0];
z q43[0];
sdg q43[0];
h q43[0];
h q43[1];
s q43[1];
z q43[1];
s q43[1];
h q43[1];
s q43[1];
sdg q43[1];
h q43[1];
y q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
z q43[0];
s q43[0];
h q43[0];
s q43[0];
h q43[1];
s q43[1];
y q43[1];
s q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
z q43[0];
h q43[0];
s q43[0];
h q43[1];
s q43[1];
z q43[1];
s q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
z q43[0];
h q43[0];
h q43[1];
s q43[1];
z q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
z q43[0];
h q43[0];
x q43[1];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
h q43[1];
s q43[1];
y q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
y q43[0];
s q43[0];
h q43[0];
s q43[0];
h q43[1];
s q43[1];
y q43[1];
h q43[1];
s q43[1];
cx q43[0],q43[1];
x q43[0];
s q43[0];
sdg q43[1];
h q43[1];
h q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
y q43[0];
h q43[0];
cx q43[0],q43[1];
h q43[0];
s q43[0];
y q43[0];
s q43[0];
h q43[1];
s q43[1];
z q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
z q43[0];
h q43[0];
h q43[1];
s q43[1];
s q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
h q43[1];
s q43[1];
cx q43[0],q43[1];
x q43[0];
z q43[1];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
x q43[0];
s q43[0];
h q43[0];
s q43[0];
sdg q43[1];
h q43[1];
y q43[1];
s q43[1];
h q43[1];
s q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
h q43[1];
s q43[1];
cx q43[0],q43[1];
x q43[0];
h q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
h q43[1];
s q43[1];
cx q43[0],q43[1];
z q43[0];
x q43[1];
s q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
z q43[0];
h q43[1];
s q43[1];
h q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
x q43[0];
s q43[0];
h q43[0];
s q43[0];
h q43[1];
s q43[1];
y q43[1];
h q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
x q43[0];
sdg q43[0];
h q43[0];
x q43[1];
h q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
h q43[0];
h q43[1];
s q43[1];
y q43[1];
h q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
z q43[0];
h q43[0];
s q43[0];
h q43[1];
s q43[1];
x q43[1];
cx q43[0],q43[1];
x q43[0];
y q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
x q43[0];
s q43[0];
h q43[0];
s q43[0];
h q43[1];
s q43[1];
cx q43[0],q43[1];
y q43[0];
sdg q43[1];
h q43[1];
x q43[1];
h q43[1];
s q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
s q43[0];
h q43[0];
s q43[0];
h q43[1];
s q43[1];
y q43[1];
h q43[1];
s q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
y q43[0];
s q43[0];
x q43[1];
s q43[1];
h q43[1];
s q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
s q43[0];
y q43[1];
s q43[1];
h q43[1];
s q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
h q43[0];
h q43[1];
s q43[1];
y q43[1];
s q43[1];
h q43[1];
s q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
z q43[0];
s q43[0];
h q43[1];
s q43[1];
y q43[1];
s q43[1];
h q43[1];
s q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
x q43[0];
sdg q43[0];
h q43[0];
sdg q43[1];
h q43[1];
z q43[1];
s q43[1];
h q43[1];
s q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
s q43[0];
h q43[0];
s q43[0];
h q43[1];
s q43[1];
z q43[1];
s q43[1];
cx q43[0],q43[1];
y q43[0];
h q43[0];
y q43[1];
h q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
x q43[0];
h q43[0];
cx q43[0],q43[1];
h q43[0];
s q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
y q43[0];
h q43[0];
s q43[0];
h q43[1];
s q43[1];
y q43[1];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
h q43[1];
s q43[1];
cx q43[0],q43[1];
y q43[0];
h q43[1];
s q43[1];
cx q43[0],q43[1];
s q43[1];
h q43[1];
s q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
x q43[0];
s q43[0];
h q43[1];
s q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
h q43[1];
s q43[1];
cx q43[0],q43[1];
h q43[0];
x q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
h q43[1];
s q43[1];
cx q43[0],q43[1];
z q43[0];
h q43[0];
s q43[0];
x q43[1];
h q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
h q43[1];
s q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
h q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
h q43[1];
s q43[1];
cx q43[0],q43[1];
y q43[0];
s q43[0];
h q43[0];
s q43[0];
x q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
y q43[0];
s q43[0];
h q43[0];
s q43[0];
y q43[0];
h q43[1];
s q43[1];
z q43[1];
sdg q43[1];
h q43[1];
z q43[1];
s q43[1];
h q43[1];
s q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
sdg q43[0];
h q43[0];
sdg q43[1];
h q43[1];
z q43[1];
s q43[1];
h q43[1];
s q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
z q43[0];
sdg q43[0];
h q43[0];
z q43[1];
s q43[1];
h q43[1];
s q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
x q43[0];
sdg q43[0];
h q43[0];
z q43[1];
h q43[1];
s q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
s q43[0];
h q43[0];
s q43[0];
sdg q43[1];
h q43[1];
x q43[1];
s q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
x q43[0];
h q43[1];
s q43[1];
z q43[1];
s q43[1];
h q43[1];
s q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
h q43[0];
s q43[0];
z q43[1];
h q43[1];
s q43[1];
cx q43[0],q43[1];
y q43[0];
h q43[0];
s q43[0];
h q43[1];
s q43[1];
s q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
h q43[1];
s q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
z q43[0];
h q43[0];
z q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
x q43[0];
x q43[1];
s q43[1];
h q43[1];
s q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
s q43[0];
h q43[1];
s q43[1];
z q43[1];
h q43[1];
s q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
h q43[1];
s q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
y q43[0];
sdg q43[0];
h q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
z q43[0];
y q43[1];
h q43[1];
s q43[1];
cx q43[0],q43[1];
y q43[0];
y q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
z q43[0];
h q43[0];
s q43[0];
h q43[1];
s q43[1];
z q43[1];
cx q43[0],q43[1];
s q43[0];
h q43[0];
s q43[0];
h q43[1];
s q43[1];
x q43[1];
s q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
h q43[1];
s q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
y q43[0];
s q43[0];
h q43[0];
s q43[0];
y q43[1];
h q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
h q43[1];
s q43[1];
cx q43[0],q43[1];
x q43[1];
h q43[1];
s q43[1];
cx q43[0],q43[1];
x q43[0];
sdg q43[0];
h q43[0];
h q43[1];
s q43[1];
y q43[1];
h q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
h q43[1];
s q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
x q43[0];
sdg q43[0];
h q43[0];
x q43[1];
s q43[1];
h q43[1];
s q43[1];
cx q43[0],q43[1];
s q43[0];
y q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
y q43[0];
sdg q43[0];
h q43[0];
h q43[1];
s q43[1];
y q43[1];
s q43[1];
h q43[1];
s q43[1];
cx q43[0],q43[1];
s q43[0];
h q43[0];
s q43[0];
x q43[1];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
x q43[0];
h q43[0];
h q43[1];
s q43[1];
z q43[1];
s q43[1];
h q43[1];
s q43[1];
cx q43[0],q43[1];
z q43[0];
h q43[1];
s q43[1];
z q43[1];
h q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
h q43[0];
sdg q43[1];
h q43[1];
y q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
x q43[0];
h q43[0];
s q43[1];
h q43[1];
s q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
y q43[0];
s q43[0];
h q43[1];
s q43[1];
z q43[1];
h q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
y q43[0];
s q43[0];
sdg q43[1];
h q43[1];
z q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
z q43[0];
h q43[0];
y q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
sdg q43[1];
h q43[1];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
h q43[1];
s q43[1];
cx q43[0],q43[1];
z q43[0];
h q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
h q43[1];
s q43[1];
cx q43[0],q43[1];
z q43[0];
h q43[0];
s q43[0];
x q43[1];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
y q43[0];
h q43[0];
sdg q43[1];
h q43[1];
x q43[1];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
z q43[0];
s q43[0];
h q43[0];
s q43[0];
h q43[1];
s q43[1];
h q43[1];
s q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
h q43[1];
s q43[1];
cx q43[0],q43[1];
z q43[0];
h q43[0];
s q43[0];
y q43[1];
h q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
h q43[1];
s q43[1];
cx q43[0],q43[1];
z q43[0];
h q43[0];
z q43[1];
s q43[1];
h q43[1];
s q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
h q43[1];
s q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
x q43[1];
h q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
h q43[1];
s q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
y q43[0];
h q43[0];
z q43[1];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
h q43[1];
s q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
y q43[1];
s q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
y q43[0];
sdg q43[0];
h q43[0];
h q43[1];
s q43[1];
y q43[1];
s q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
h q43[1];
s q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
z q43[0];
h q43[0];
z q43[1];
h q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
x q43[0];
s q43[0];
h q43[0];
s q43[0];
y q43[1];
h q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
x q43[0];
s q43[0];
h q43[0];
s q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
s q43[0];
h q43[0];
s q43[0];
h q43[1];
s q43[1];
z q43[1];
h q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
y q43[0];
h q43[1];
s q43[1];
z q43[1];
s q43[1];
cx q43[0],q43[1];
x q43[0];
h q43[1];
s q43[1];
x q43[1];
s q43[1];
h q43[1];
s q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
y q43[0];
s q43[0];
h q43[1];
s q43[1];
s q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
z q43[0];
sdg q43[0];
h q43[0];
h q43[1];
s q43[1];
y q43[1];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
y q43[0];
h q43[0];
s q43[0];
h q43[1];
s q43[1];
x q43[1];
h q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
x q43[0];
s q43[0];
h q43[1];
s q43[1];
y q43[1];
h q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
x q43[0];
s q43[0];
h q43[0];
s q43[0];
x q43[1];
s q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
x q43[0];
s q43[0];
h q43[1];
s q43[1];
s q43[1];
h q43[1];
s q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
z q43[0];
z q43[1];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
z q43[0];
sdg q43[0];
h q43[0];
z q43[1];
s q43[1];
h q43[1];
s q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
z q43[0];
h q43[0];
x q43[1];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
z q43[0];
h q43[0];
s q43[0];
sdg q43[1];
h q43[1];
y q43[1];
h q43[1];
s q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
z q43[0];
h q43[0];
h q43[1];
s q43[1];
h q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
y q43[0];
s q43[0];
h q43[0];
s q43[0];
h q43[1];
s q43[1];
x q43[1];
s q43[1];
h q43[1];
s q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
y q43[0];
h q43[0];
h q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
x q43[0];
sdg q43[0];
h q43[0];
z q43[1];
s q43[1];
cx q43[0],q43[1];
x q43[0];
h q43[0];
s q43[0];
x q43[1];
s q43[1];
h q43[1];
s q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
x q43[0];
h q43[1];
s q43[1];
x q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
y q43[0];
sdg q43[0];
h q43[0];
sdg q43[1];
h q43[1];
x q43[1];
h q43[1];
cx q43[0],q43[1];
z q43[0];
s q43[0];
h q43[1];
s q43[1];
y q43[1];
s q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
x q43[0];
h q43[0];
s q43[0];
h q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
z q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
x q43[0];
h q43[0];
h q43[1];
s q43[1];
y q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
y q43[0];
sdg q43[0];
h q43[0];
h q43[1];
s q43[1];
z q43[1];
h q43[1];
s q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
h q43[1];
s q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
sdg q43[0];
h q43[0];
z q43[1];
s q43[1];
h q43[1];
s q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
z q43[0];
sdg q43[0];
h q43[0];
x q43[1];
h q43[1];
s q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
h q43[1];
s q43[1];
cx q43[0],q43[1];
x q43[0];
x q43[1];
h q43[1];
cx q43[0],q43[1];
x q43[1];
h q43[1];
s q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
h q43[1];
s q43[1];
cx q43[0],q43[1];
y q43[0];
y q43[1];
h q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
h q43[1];
s q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
sdg q43[0];
h q43[0];
z q43[1];
h q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
s q43[0];
h q43[0];
s q43[0];
z q43[1];
s q43[1];
h q43[1];
s q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
h q43[1];
s q43[1];
cx q43[0],q43[1];
y q43[0];
sdg q43[0];
h q43[0];
s q43[1];
h q43[1];
s q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
sdg q43[0];
h q43[0];
sdg q43[1];
h q43[1];
z q43[1];
s q43[1];
cx q43[0],q43[1];
s q43[0];
h q43[0];
s q43[0];
sdg q43[0];
h q43[0];
h q43[1];
s q43[1];
x q43[1];
h q43[1];
s q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
z q43[0];
h q43[0];
s q43[0];
x q43[1];
h q43[1];
s q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
y q43[1];
s q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
h q43[1];
s q43[1];
cx q43[0],q43[1];
z q43[0];
s q43[0];
x q43[1];
cx q43[0],q43[1];
z q43[0];
sdg q43[1];
h q43[1];
h q43[1];
s q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
x q43[0];
sdg q43[0];
h q43[0];
sdg q43[1];
h q43[1];
z q43[1];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
y q43[0];
h q43[0];
s q43[0];
cx q43[0],q43[1];
h q43[0];
s q43[0];
s q43[0];
h q43[0];
s q43[0];
h q43[1];
s q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
y q43[0];
h q43[0];
s q43[0];
h q43[1];
s q43[1];
x q43[1];
h q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
z q43[0];
h q43[1];
s q43[1];
y q43[1];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
z q43[0];
h q43[1];
s q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
x q43[0];
h q43[0];
h q43[1];
s q43[1];
s q43[1];
h q43[1];
s q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
h q43[1];
s q43[1];
cx q43[0],q43[1];
x q43[0];
sdg q43[0];
h q43[0];
x q43[1];
h q43[1];
s q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
h q43[1];
s q43[1];
cx q43[0],q43[1];
y q43[0];
s q43[0];
h q43[0];
s q43[0];
s q43[1];
h q43[1];
s q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
z q43[0];
h q43[0];
h q43[1];
s q43[1];
x q43[1];
s q43[1];
h q43[1];
s q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
z q43[0];
h q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
z q43[0];
sdg q43[1];
h q43[1];
x q43[1];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
z q43[0];
h q43[0];
sdg q43[1];
h q43[1];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
y q43[0];
h q43[0];
s q43[0];
y q43[1];
s q43[1];
h q43[1];
s q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
h q43[1];
s q43[1];
x q43[1];
s q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
z q43[0];
s q43[0];
h q43[0];
s q43[0];
z q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
x q43[0];
h q43[0];
y q43[1];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
z q43[0];
h q43[0];
y q43[1];
s q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
y q43[0];
h q43[0];
sdg q43[1];
h q43[1];
y q43[1];
h q43[1];
s q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
h q43[1];
s q43[1];
y q43[1];
s q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
h q43[1];
s q43[1];
cx q43[0],q43[1];
x q43[0];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
z q43[0];
sdg q43[0];
h q43[0];
h q43[1];
s q43[1];
h q43[1];
x q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
s q43[0];
h q43[0];
s q43[0];
y q43[1];
h q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
h q43[1];
s q43[1];
cx q43[0],q43[1];
s q43[0];
h q43[0];
s q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
s q43[0];
h q43[1];
s q43[1];
z q43[1];
s q43[1];
h q43[1];
s q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
z q43[0];
sdg q43[0];
h q43[0];
z q43[1];
h q43[1];
s q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
h q43[1];
s q43[1];
cx q43[0],q43[1];
s q43[0];
x q43[1];
h q43[1];
s q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
sdg q43[0];
h q43[0];
h q43[1];
s q43[1];
y q43[1];
s q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
x q43[0];
s q43[0];
x q43[1];
h q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
x q43[0];
s q43[0];
y q43[1];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
y q43[0];
h q43[0];
x q43[1];
s q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
h q43[1];
s q43[1];
x q43[1];
s q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
y q43[0];
s q43[0];
h q43[0];
s q43[0];
h q43[1];
s q43[1];
x q43[1];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
z q43[0];
h q43[1];
s q43[1];
y q43[1];
h q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
x q43[0];
h q43[1];
s q43[1];
x q43[1];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
y q43[0];
z q43[1];
s q43[1];
h q43[1];
s q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
y q43[0];
sdg q43[0];
h q43[0];
h q43[1];
s q43[1];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
y q43[1];
h q43[1];
s q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
z q43[0];
h q43[0];
h q43[1];
s q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
x q43[0];
sdg q43[0];
h q43[0];
y q43[0];
s q43[0];
h q43[0];
s q43[0];
h q43[1];
s q43[1];
y q43[1];
x q43[1];
h q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
z q43[0];
y q43[1];
h q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
x q43[0];
sdg q43[1];
h q43[1];
x q43[1];
s q43[1];
h q43[1];
s q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
z q43[0];
sdg q43[0];
h q43[0];
h q43[1];
cx q43[0],q43[1];
x q43[0];
sdg q43[1];
h q43[1];
y q43[1];
h q43[1];
s q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
y q43[0];
sdg q43[1];
h q43[1];
y q43[1];
h q43[1];
cx q43[0],q43[1];
s q43[0];
h q43[0];
s q43[0];
h q43[1];
s q43[1];
x q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
h q43[1];
s q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
y q43[0];
h q43[0];
s q43[0];
x q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
s q43[0];
z q43[1];
s q43[1];
h q43[1];
s q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
h q43[1];
s q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
y q43[0];
s q43[0];
h q43[0];
s q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
h q43[1];
s q43[1];
cx q43[0],q43[1];
h q43[0];
h q43[1];
s q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
y q43[0];
z q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
x q43[0];
h q43[0];
s q43[0];
sdg q43[1];
h q43[1];
y q43[1];
h q43[1];
s q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
x q43[0];
s q43[0];
h q43[1];
s q43[1];
s q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
h q43[1];
s q43[1];
cx q43[0],q43[1];
h q43[0];
h q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
h q43[0];
s q43[0];
h q43[1];
s q43[1];
y q43[1];
s q43[1];
cx q43[0],q43[1];
h q43[0];
x q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
x q43[0];
sdg q43[0];
h q43[0];
sdg q43[1];
h q43[1];
y q43[1];
s q43[1];
h q43[1];
s q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
y q43[0];
h q43[0];
s q43[0];
h q43[1];
s q43[1];
z q43[1];
s q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
z q43[0];
h q43[0];
h q43[1];
s q43[1];
h q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
x q43[0];
s q43[0];
h q43[0];
s q43[0];
h q43[1];
s q43[1];
y q43[1];
s q43[1];
h q43[1];
s q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
h q43[1];
s q43[1];
cx q43[0],q43[1];
y q43[0];
h q43[0];
s q43[0];
y q43[1];
h q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
y q43[0];
z q43[1];
h q43[1];
s q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
y q43[0];
h q43[1];
s q43[1];
y q43[1];
s q43[1];
h q43[1];
s q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
h q43[0];
s q43[0];
h q43[1];
s q43[1];
s q43[1];
h q43[1];
s q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
x q43[0];
s q43[0];
y q43[1];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
x q43[0];
h q43[0];
s q43[0];
h q43[1];
s q43[1];
y q43[1];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
h q43[1];
s q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
z q43[1];
h q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
x q43[0];
x q43[1];
s q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
h q43[0];
h q43[1];
s q43[1];
x q43[1];
h q43[1];
cx q43[0],q43[1];
y q43[0];
h q43[0];
sdg q43[1];
h q43[1];
x q43[1];
s q43[1];
h q43[1];
s q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
x q43[0];
h q43[0];
s q43[0];
h q43[1];
s q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
h q43[1];
s q43[1];
cx q43[0],q43[1];
y q43[0];
sdg q43[0];
h q43[0];
h q43[1];
s q43[1];
cx q43[0],q43[1];
z q43[0];
h q43[0];
h q43[1];
s q43[1];
z q43[1];
s q43[1];
cx q43[0],q43[1];
y q43[0];
h q43[0];
h q43[1];
s q43[1];
z q43[1];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
z q43[0];
s q43[0];
sdg q43[1];
h q43[1];
y q43[1];
h q43[1];
cx q43[0],q43[1];
y q43[0];
h q43[0];
y q43[1];
s q43[1];
h q43[1];
s q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
h q43[0];
s q43[0];
h q43[1];
s q43[1];
z q43[1];
s q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
x q43[0];
s q43[0];
h q43[0];
s q43[0];
z q43[0];
h q43[0];
y q43[1];
s q43[1];
h q43[1];
s q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
z q43[0];
h q43[0];
sdg q43[1];
h q43[1];
y q43[1];
s q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
z q43[0];
h q43[0];
h q43[1];
s q43[1];
y q43[1];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
h q43[1];
s q43[1];
cx q43[0],q43[1];
x q43[0];
s q43[0];
h q43[0];
s q43[0];
x q43[1];
s q43[1];
h q43[1];
s q43[1];
cx q43[0],q43[1];
z q43[0];
h q43[0];
s q43[0];
y q43[1];
s q43[1];
h q43[1];
s q43[1];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
x q43[0];
s q43[0];
h q43[0];
s q43[0];
z q43[1];
h q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
h q43[0];
s q43[0];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
x q43[0];
sdg q43[0];
h q43[0];
z q43[1];
h q43[1];
s q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
y q43[0];
s q43[0];
x q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
h q43[1];
s q43[1];
cx q43[0],q43[1];
z q43[0];
s q43[0];
h q43[0];
s q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
h q43[1];
s q43[1];
cx q43[0],q43[1];
z q43[0];
s q43[0];
z q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
z q43[0];
h q43[0];
x q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
x q43[0];
s q43[0];
h q43[0];
s q43[0];
h q43[1];
s q43[1];
y q43[1];
s q43[1];
h q43[1];
s q43[1];
cx q43[0],q43[1];
h q43[0];
y q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
x q43[0];
h q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
h q43[1];
s q43[1];
cx q43[0],q43[1];
x q43[0];
sdg q43[0];
h q43[0];
x q43[1];
s q43[1];
h q43[1];
s q43[1];
cx q43[0],q43[1];
z q43[0];
sdg q43[0];
h q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
h q43[0];
h q43[1];
s q43[1];
y q43[1];
h q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
y q43[0];
s q43[0];
h q43[0];
s q43[0];
sdg q43[1];
h q43[1];
y q43[1];
h q43[1];
s q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
y q43[0];
sdg q43[1];
h q43[1];
h q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
y q43[0];
h q43[0];
s q43[0];
s q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
s q43[0];
s q43[1];
h q43[1];
s q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
z q43[0];
h q43[0];
s q43[0];
sdg q43[1];
h q43[1];
s q43[1];
h q43[1];
s q43[1];
cx q43[0],q43[1];
h q43[0];
h q43[1];
s q43[1];
h q43[1];
s q43[1];
cx q43[0],q43[1];
x q43[0];
h q43[1];
s q43[1];
h q43[1];
s q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
y q43[0];
sdg q43[0];
h q43[0];
s q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
h q43[1];
s q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
z q43[0];
x q43[1];
s q43[1];
h q43[1];
s q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
x q43[0];
h q43[0];
h q43[1];
s q43[1];
y q43[1];
h q43[1];
s q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
x q43[0];
h q43[0];
s q43[0];
h q43[0];
h q43[1];
s q43[1];
z q43[1];
h q43[1];
s q43[1];
h q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
z q43[0];
h q43[0];
s q43[0];
y q43[1];
s q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
y q43[0];
s q43[0];
x q43[1];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
x q43[0];
sdg q43[0];
h q43[0];
y q43[1];
h q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
x q43[0];
h q43[0];
s q43[0];
sdg q43[1];
h q43[1];
z q43[1];
s q43[1];
h q43[1];
s q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
z q43[0];
h q43[1];
s q43[1];
h q43[1];
s q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
h q43[0];
s q43[0];
h q43[1];
s q43[1];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
y q43[0];
h q43[0];
h q43[1];
s q43[1];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
y q43[0];
s q43[0];
h q43[0];
s q43[0];
h q43[1];
s q43[1];
y q43[1];
s q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
h q43[1];
s q43[1];
cx q43[0],q43[1];
z q43[0];
sdg q43[0];
h q43[0];
x q43[1];
s q43[1];
h q43[1];
s q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
x q43[0];
s q43[0];
h q43[1];
s q43[1];
h q43[1];
s q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
x q43[0];
sdg q43[1];
h q43[1];
x q43[1];
h q43[1];
cx q43[0],q43[1];
z q43[0];
sdg q43[1];
h q43[1];
s q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
h q43[1];
s q43[1];
cx q43[0],q43[1];
y q43[0];
h q43[0];
y q43[1];
s q43[1];
h q43[1];
s q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
h q43[1];
s q43[1];
cx q43[0],q43[1];
y q43[0];
h q43[0];
z q43[1];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
h q43[1];
s q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
s q43[0];
y q43[1];
s q43[1];
h q43[1];
s q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
h q43[0];
x q43[0];
h q43[0];
s q43[0];
x q43[1];
h q43[1];
h q43[1];
s q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
h q43[1];
s q43[1];
cx q43[0],q43[1];
z q43[0];
s q43[0];
h q43[0];
s q43[0];
z q43[0];
h q43[0];
z q43[1];
h q43[1];
z q43[1];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
sdg q43[0];
h q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
h q43[1];
s q43[1];
cx q43[0],q43[1];
z q43[0];
s q43[0];
h q43[0];
s q43[0];
y q43[1];
s q43[1];
h q43[1];
s q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
h q43[1];
s q43[1];
cx q43[0],q43[1];
z q43[0];
y q43[1];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
z q43[0];
h q43[0];
s q43[0];
h q43[1];
s q43[1];
s q43[1];
h q43[1];
s q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
z q43[0];
h q43[1];
s q43[1];
z q43[1];
h q43[1];
s q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
h q43[1];
s q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
y q43[0];
sdg q43[0];
h q43[0];
s q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
x q43[0];
h q43[1];
s q43[1];
z q43[1];
s q43[1];
h q43[1];
s q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
y q43[0];
s q43[0];
h q43[0];
s q43[0];
h q43[1];
s q43[1];
y q43[1];
s q43[1];
h q43[1];
s q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
y q43[0];
h q43[0];
h q43[1];
s q43[1];
y q43[1];
s q43[1];
h q43[1];
s q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
z q43[1];
h q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
h q43[1];
s q43[1];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
x q43[0];
h q43[1];
s q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
x q43[0];
h q43[0];
s q43[0];
y q43[1];
s q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
h q43[1];
s q43[1];
cx q43[0],q43[1];
s q43[0];
x q43[1];
h q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
h q43[1];
s q43[1];
cx q43[0],q43[1];
x q43[0];
y q43[1];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
z q43[0];
sdg q43[0];
h q43[0];
h q43[1];
s q43[1];
z q43[1];
s q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
x q43[0];
h q43[0];
s q43[0];
h q43[1];
s q43[1];
z q43[1];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
x q43[0];
y q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
x q43[0];
sdg q43[0];
h q43[0];
sdg q43[1];
h q43[1];
y q43[1];
h q43[1];
s q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
x q43[0];
sdg q43[0];
h q43[0];
h q43[1];
s q43[1];
z q43[1];
s q43[1];
h q43[1];
s q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
y q43[1];
s q43[1];
h q43[1];
s q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
h q43[1];
s q43[1];
cx q43[0],q43[1];
z q43[0];
h q43[0];
s q43[0];
z q43[1];
h q43[1];
cx q43[0],q43[1];
z q43[0];
s q43[0];
y q43[1];
h q43[1];
s q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
z q43[0];
h q43[0];
s q43[0];
x q43[1];
s q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
z q43[0];
h q43[0];
z q43[1];
cx q43[0],q43[1];
x q43[0];
h q43[0];
h q43[1];
s q43[1];
x q43[1];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
z q43[0];
sdg q43[0];
h q43[0];
y q43[1];
s q43[1];
h q43[1];
s q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
s q43[0];
h q43[0];
s q43[0];
y q43[0];
s q43[1];
y q43[1];
s q43[1];
h q43[1];
s q43[1];
cx q43[0],q43[1];
z q43[0];
h q43[1];
s q43[1];
x q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
x q43[0];
s q43[0];
h q43[0];
s q43[0];
z q43[1];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
h q43[0];
s q43[0];
z q43[1];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
h q43[0];
h q43[1];
s q43[1];
y q43[1];
h q43[1];
s q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
z q43[0];
sdg q43[0];
h q43[0];
y q43[1];
cx q43[0],q43[1];
y q43[0];
s q43[0];
h q43[1];
s q43[1];
y q43[1];
h q43[1];
s q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
x q43[0];
sdg q43[0];
h q43[0];
z q43[0];
h q43[0];
s q43[0];
h q43[1];
s q43[1];
y q43[1];
y q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
h q43[1];
s q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
s q43[0];
h q43[0];
s q43[0];
y q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
z q43[0];
s q43[0];
h q43[0];
s q43[0];
z q43[1];
s q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
h q43[1];
s q43[1];
cx q43[0],q43[1];
z q43[0];
cx q43[0],q43[1];
h q43[1];
s q43[1];
x q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
s q43[0];
h q43[0];
s q43[0];
sdg q43[1];
h q43[1];
y q43[1];
s q43[1];
h q43[1];
s q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
s q43[0];
x q43[1];
s q43[1];
h q43[1];
s q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
h q43[1];
s q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
z q43[0];
h q43[0];
s q43[0];
s q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
z q43[0];
h q43[0];
s q43[0];
h q43[1];
s q43[1];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
h q43[1];
s q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
z q43[0];
h q43[0];
s q43[0];
x q43[1];
s q43[1];
h q43[1];
s q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
h q43[1];
s q43[1];
cx q43[0],q43[1];
y q43[0];
s q43[0];
z q43[1];
h q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
h q43[1];
s q43[1];
cx q43[0],q43[1];
s q43[0];
h q43[0];
s q43[0];
x q43[1];
s q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
s q43[0];
h q43[1];
s q43[1];
h q43[1];
s q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
x q43[0];
h q43[0];
h q43[1];
s q43[1];
z q43[1];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
y q43[0];
h q43[0];
s q43[0];
h q43[1];
s q43[1];
z q43[1];
s q43[1];
y q43[1];
s q43[1];
h q43[1];
s q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
y q43[0];
s q43[0];
z q43[1];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
h q43[1];
s q43[1];
cx q43[0],q43[1];
x q43[0];
s q43[0];
h q43[0];
s q43[0];
z q43[1];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
z q43[0];
y q43[1];
h q43[1];
cx q43[0],q43[1];
x q43[0];
h q43[0];
s q43[0];
y q43[1];
s q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
z q43[0];
h q43[0];
s q43[0];
y q43[1];
h q43[1];
s q43[1];
cx q43[0],q43[1];
z q43[0];
sdg q43[0];
h q43[0];
sdg q43[1];
h q43[1];
y q43[1];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
h q43[0];
s q43[0];
h q43[1];
s q43[1];
cx q43[0],q43[1];
y q43[0];
h q43[0];
s q43[0];
y q43[1];
sdg q43[1];
h q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
h q43[1];
s q43[1];
cx q43[0],q43[1];
s q43[0];
x q43[1];
s q43[1];
cx q43[0],q43[1];
sdg q43[0];
h q43[0];
sdg q43[1];
h q43[1];
