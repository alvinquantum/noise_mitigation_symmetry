OPENQASM 2.0;
include "qelib1.inc";
qreg q35[2];
y q35[0];
h q35[0];
h q35[1];
s q35[1];
cx q35[0],q35[1];
h q35[0];
s q35[0];
x q35[0];
s q35[0];
h q35[1];
s q35[1];
y q35[1];
h q35[1];
s q35[1];
cx q35[0],q35[1];
sdg q35[0];
h q35[0];
sdg q35[1];
h q35[1];
cx q35[0],q35[1];
z q35[0];
h q35[0];
y q35[1];
cx q35[0],q35[1];
sdg q35[0];
h q35[0];
h q35[1];
s q35[1];
cx q35[0],q35[1];
x q35[0];
sdg q35[0];
h q35[0];
x q35[0];
h q35[0];
s q35[0];
z q35[1];
s q35[1];
cx q35[0],q35[1];
sdg q35[0];
h q35[0];
h q35[1];
s q35[1];
cx q35[0],q35[1];
h q35[0];
s q35[0];
x q35[0];
s q35[0];
h q35[0];
s q35[0];
y q35[1];
cx q35[0],q35[1];
h q35[0];
s q35[0];
x q35[0];
h q35[1];
s q35[1];
cx q35[0],q35[1];
h q35[0];
s q35[0];
y q35[0];
s q35[0];
y q35[1];
sdg q35[1];
h q35[1];
cx q35[0],q35[1];
h q35[0];
s q35[0];
sdg q35[1];
h q35[1];
cx q35[0],q35[1];
s q35[0];
h q35[0];
s q35[0];
h q35[1];
s q35[1];
y q35[1];
sdg q35[1];
h q35[1];
cx q35[0],q35[1];
sdg q35[0];
h q35[0];
z q35[0];
s q35[0];
y q35[1];
cx q35[0],q35[1];
sdg q35[0];
h q35[0];
sdg q35[1];
h q35[1];
cx q35[0],q35[1];
x q35[0];
h q35[0];
cx q35[0],q35[1];
sdg q35[0];
h q35[0];
sdg q35[1];
h q35[1];
cx q35[0],q35[1];
h q35[0];
s q35[0];
h q35[1];
s q35[1];
cx q35[0],q35[1];
h q35[0];
s q35[0];
sdg q35[1];
h q35[1];
cx q35[0],q35[1];
sdg q35[0];
h q35[0];
s q35[0];
h q35[0];
s q35[0];
sdg q35[1];
h q35[1];
y q35[1];
s q35[1];
cx q35[0],q35[1];
h q35[0];
s q35[0];
y q35[0];
x q35[1];
cx q35[0],q35[1];
sdg q35[0];
h q35[0];
h q35[1];
s q35[1];
cx q35[0],q35[1];
h q35[0];
s q35[0];
z q35[0];
s q35[0];
s q35[1];
cx q35[0],q35[1];
h q35[0];
s q35[0];
x q35[0];
h q35[0];
s q35[0];
h q35[1];
s q35[1];
y q35[1];
h q35[1];
cx q35[0],q35[1];
sdg q35[0];
h q35[0];
x q35[0];
sdg q35[0];
h q35[0];
sdg q35[1];
h q35[1];
h q35[1];
s q35[1];
cx q35[0],q35[1];
sdg q35[0];
h q35[0];
sdg q35[1];
h q35[1];
cx q35[0],q35[1];
z q35[1];
cx q35[0],q35[1];
sdg q35[0];
h q35[0];
sdg q35[1];
h q35[1];
cx q35[0],q35[1];
y q35[0];
h q35[0];
h q35[1];
s q35[1];
x q35[1];
s q35[1];
cx q35[0],q35[1];
sdg q35[0];
h q35[0];
sdg q35[1];
h q35[1];
cx q35[0],q35[1];
x q35[0];
x q35[1];
h q35[1];
s q35[1];
cx q35[0],q35[1];
sdg q35[0];
h q35[0];
x q35[0];
sdg q35[0];
h q35[0];
sdg q35[1];
h q35[1];
h q35[1];
cx q35[0],q35[1];
h q35[0];
s q35[0];
h q35[1];
s q35[1];
cx q35[0],q35[1];
x q35[0];
h q35[0];
s q35[0];
y q35[1];
h q35[1];
cx q35[0],q35[1];
x q35[0];
s q35[0];
h q35[0];
s q35[0];
sdg q35[1];
h q35[1];
x q35[1];
h q35[1];
s q35[1];
cx q35[0],q35[1];
h q35[0];
s q35[0];
sdg q35[1];
h q35[1];
cx q35[0],q35[1];
x q35[0];
h q35[0];
h q35[1];
s q35[1];
cx q35[0],q35[1];
sdg q35[0];
h q35[0];
h q35[1];
s q35[1];
cx q35[0],q35[1];
y q35[1];
s q35[1];
h q35[1];
s q35[1];
cx q35[0],q35[1];
sdg q35[0];
h q35[0];
x q35[0];
h q35[0];
s q35[0];
h q35[1];
s q35[1];
s q35[1];
h q35[1];
s q35[1];
cx q35[0],q35[1];
sdg q35[0];
h q35[0];
sdg q35[1];
h q35[1];
cx q35[0],q35[1];
h q35[1];
s q35[1];
z q35[1];
h q35[1];
s q35[1];
cx q35[0],q35[1];
x q35[0];
sdg q35[0];
h q35[0];
h q35[1];
s q35[1];
z q35[1];
s q35[1];
h q35[1];
s q35[1];
cx q35[0],q35[1];
h q35[0];
s q35[0];
y q35[0];
sdg q35[0];
h q35[0];
s q35[0];
y q35[1];
h q35[1];
s q35[1];
s q35[1];
h q35[1];
s q35[1];
cx q35[0],q35[1];
x q35[0];
sdg q35[0];
h q35[0];
h q35[1];
s q35[1];
sdg q35[1];
h q35[1];
cx q35[0],q35[1];
sdg q35[0];
h q35[0];
sdg q35[1];
h q35[1];
cx q35[0],q35[1];
h q35[0];
s q35[0];
y q35[0];
sdg q35[0];
h q35[0];
x q35[1];
h q35[1];
cx q35[0],q35[1];
h q35[0];
y q35[1];
h q35[1];
s q35[1];
cx q35[0],q35[1];
h q35[0];
s q35[0];
sdg q35[1];
h q35[1];
cx q35[0],q35[1];
s q35[0];
h q35[0];
s q35[0];
z q35[1];
sdg q35[1];
h q35[1];
cx q35[0],q35[1];
h q35[0];
s q35[0];
sdg q35[1];
h q35[1];
cx q35[0],q35[1];
z q35[0];
sdg q35[0];
h q35[0];
s q35[1];
cx q35[0],q35[1];
x q35[0];
h q35[0];
h q35[1];
s q35[1];
s q35[1];
cx q35[0],q35[1];
y q35[0];
sdg q35[0];
h q35[0];
h q35[1];
s q35[1];
z q35[1];
sdg q35[1];
h q35[1];
cx q35[0],q35[1];
sdg q35[0];
h q35[0];
sdg q35[1];
h q35[1];
cx q35[0],q35[1];
h q35[0];
s q35[0];
x q35[0];
s q35[0];
h q35[1];
s q35[1];
z q35[1];
h q35[1];
s q35[1];
cx q35[0],q35[1];
sdg q35[0];
h q35[0];
sdg q35[1];
h q35[1];
cx q35[0],q35[1];
h q35[0];
s q35[0];
h q35[1];
s q35[1];
cx q35[0],q35[1];
h q35[0];
s q35[0];
z q35[1];
h q35[1];
cx q35[0],q35[1];
sdg q35[0];
h q35[0];
x q35[0];
s q35[0];
sdg q35[1];
h q35[1];
z q35[1];
h q35[1];
cx q35[0],q35[1];
sdg q35[0];
h q35[0];
sdg q35[1];
h q35[1];
cx q35[0],q35[1];
s q35[0];
h q35[0];
s q35[0];
x q35[1];
sdg q35[1];
h q35[1];
cx q35[0],q35[1];
h q35[0];
s q35[0];
y q35[0];
h q35[0];
x q35[1];
h q35[1];
s q35[1];
cx q35[0],q35[1];
sdg q35[0];
h q35[0];
sdg q35[1];
h q35[1];
cx q35[0],q35[1];
h q35[0];
s q35[0];
z q35[0];
s q35[0];
h q35[0];
s q35[0];
y q35[1];
s q35[1];
cx q35[0],q35[1];
sdg q35[0];
h q35[0];
h q35[1];
s q35[1];
cx q35[0],q35[1];
h q35[0];
s q35[0];
y q35[0];
sdg q35[0];
h q35[0];
y q35[1];
cx q35[0],q35[1];
sdg q35[0];
h q35[0];
sdg q35[1];
h q35[1];
cx q35[0],q35[1];
h q35[0];
s q35[0];
x q35[0];
sdg q35[0];
h q35[0];
h q35[1];
s q35[1];
z q35[1];
h q35[1];
cx q35[0],q35[1];
sdg q35[0];
h q35[0];
sdg q35[1];
h q35[1];
cx q35[0],q35[1];
h q35[0];
s q35[0];
sdg q35[0];
h q35[0];
y q35[1];
h q35[1];
cx q35[0],q35[1];
h q35[0];
s q35[0];
y q35[0];
h q35[0];
sdg q35[1];
h q35[1];
h q35[1];
cx q35[0],q35[1];
sdg q35[0];
h q35[0];
sdg q35[1];
h q35[1];
cx q35[0],q35[1];
h q35[0];
s q35[0];
y q35[1];
h q35[1];
cx q35[0],q35[1];
h q35[0];
s q35[0];
h q35[1];
s q35[1];
cx q35[0],q35[1];
z q35[0];
sdg q35[0];
h q35[0];
x q35[1];
s q35[1];
h q35[1];
s q35[1];
cx q35[0],q35[1];
sdg q35[0];
h q35[0];
z q35[0];
h q35[0];
s q35[0];
sdg q35[1];
h q35[1];
cx q35[0],q35[1];
sdg q35[0];
h q35[0];
sdg q35[1];
h q35[1];
cx q35[0],q35[1];
y q35[0];
h q35[0];
s q35[0];
h q35[1];
s q35[1];
cx q35[0],q35[1];
sdg q35[0];
h q35[0];
h q35[1];
s q35[1];
cx q35[0],q35[1];
h q35[0];
s q35[0];
h q35[0];
s q35[0];
y q35[1];
s q35[1];
cx q35[0],q35[1];
sdg q35[0];
h q35[0];
h q35[1];
s q35[1];
cx q35[0],q35[1];
z q35[0];
h q35[0];
s q35[0];
x q35[1];
sdg q35[1];
h q35[1];
cx q35[0],q35[1];
h q35[0];
s q35[0];
sdg q35[1];
h q35[1];
cx q35[0],q35[1];
z q35[0];
sdg q35[0];
h q35[0];
h q35[1];
s q35[1];
x q35[1];
s q35[1];
cx q35[0],q35[1];
sdg q35[0];
h q35[0];
sdg q35[1];
h q35[1];
cx q35[0],q35[1];
h q35[0];
s q35[0];
h q35[1];
s q35[1];
cx q35[0],q35[1];
x q35[0];
h q35[0];
s q35[0];
x q35[1];
s q35[1];
cx q35[0],q35[1];
h q35[0];
s q35[0];
s q35[0];
h q35[0];
s q35[0];
s q35[0];
h q35[0];
s q35[0];
sdg q35[1];
h q35[1];
z q35[1];
s q35[1];
h q35[1];
s q35[1];
z q35[1];
sdg q35[1];
h q35[1];
cx q35[0],q35[1];
h q35[0];
s q35[0];
sdg q35[0];
h q35[0];
sdg q35[1];
h q35[1];
z q35[1];
h q35[1];
s q35[1];
cx q35[0],q35[1];
h q35[0];
s q35[0];
sdg q35[1];
h q35[1];
cx q35[0],q35[1];
y q35[0];
z q35[1];
cx q35[0],q35[1];
sdg q35[0];
h q35[0];
y q35[0];
s q35[0];
h q35[0];
s q35[0];
x q35[1];
sdg q35[1];
h q35[1];
cx q35[0],q35[1];
sdg q35[0];
h q35[0];
h q35[1];
s q35[1];
cx q35[0],q35[1];
x q35[0];
h q35[0];
z q35[1];
h q35[1];
s q35[1];
cx q35[0],q35[1];
sdg q35[0];
h q35[0];
h q35[1];
s q35[1];
cx q35[0],q35[1];
x q35[0];
sdg q35[0];
h q35[0];
y q35[1];
h q35[1];
s q35[1];
cx q35[0],q35[1];
sdg q35[0];
h q35[0];
sdg q35[1];
h q35[1];
cx q35[0],q35[1];
s q35[0];
h q35[1];
s q35[1];
y q35[1];
h q35[1];
cx q35[0],q35[1];
sdg q35[0];
h q35[0];
h q35[1];
s q35[1];
cx q35[0],q35[1];
x q35[0];
x q35[1];
sdg q35[1];
h q35[1];
cx q35[0],q35[1];
h q35[0];
s q35[0];
h q35[1];
s q35[1];
cx q35[0],q35[1];
y q35[0];
s q35[0];
h q35[0];
s q35[0];
x q35[1];
cx q35[0],q35[1];
sdg q35[0];
h q35[0];
z q35[0];
h q35[1];
s q35[1];
x q35[1];
sdg q35[1];
h q35[1];
cx q35[0],q35[1];
x q35[0];
s q35[0];
h q35[1];
s q35[1];
s q35[1];
cx q35[0],q35[1];
h q35[0];
s q35[0];
sdg q35[1];
h q35[1];
y q35[1];
s q35[1];
cx q35[0],q35[1];
sdg q35[0];
h q35[0];
sdg q35[1];
h q35[1];
cx q35[0],q35[1];
y q35[0];
sdg q35[0];
h q35[0];
sdg q35[1];
h q35[1];
cx q35[0],q35[1];
h q35[0];
s q35[0];
h q35[1];
s q35[1];
cx q35[0],q35[1];
z q35[0];
z q35[1];
h q35[1];
cx q35[0],q35[1];
sdg q35[0];
h q35[0];
sdg q35[1];
h q35[1];
cx q35[0],q35[1];
h q35[0];
s q35[0];
h q35[1];
s q35[1];
cx q35[0],q35[1];
h q35[0];
cx q35[0],q35[1];
sdg q35[0];
h q35[0];
y q35[0];
h q35[0];
sdg q35[1];
h q35[1];
z q35[1];
sdg q35[1];
h q35[1];
cx q35[0],q35[1];
sdg q35[0];
h q35[0];
sdg q35[1];
h q35[1];
cx q35[0],q35[1];
z q35[0];
sdg q35[0];
h q35[0];
h q35[1];
s q35[1];
x q35[1];
s q35[1];
h q35[1];
s q35[1];
cx q35[0],q35[1];
h q35[0];
s q35[0];
x q35[0];
s q35[0];
h q35[0];
s q35[0];
x q35[1];
h q35[1];
s q35[1];
cx q35[0],q35[1];
y q35[0];
s q35[0];
h q35[0];
s q35[0];
h q35[1];
s q35[1];
sdg q35[1];
h q35[1];
cx q35[0],q35[1];
h q35[0];
s q35[0];
h q35[1];
s q35[1];
cx q35[0],q35[1];
s q35[0];
h q35[0];
s q35[0];
z q35[0];
s q35[0];
h q35[0];
s q35[0];
z q35[1];
y q35[1];
s q35[1];
h q35[1];
s q35[1];
cx q35[0],q35[1];
y q35[0];
sdg q35[1];
h q35[1];
y q35[1];
cx q35[0],q35[1];
sdg q35[0];
h q35[0];
sdg q35[1];
h q35[1];
cx q35[0],q35[1];
z q35[0];
y q35[1];
cx q35[0],q35[1];
h q35[0];
s q35[0];
sdg q35[1];
h q35[1];
cx q35[0],q35[1];
y q35[0];
h q35[0];
s q35[0];
x q35[0];
h q35[0];
h q35[1];
s q35[1];
z q35[1];
h q35[1];
x q35[1];
cx q35[0],q35[1];
x q35[0];
h q35[1];
s q35[1];
y q35[1];
s q35[1];
cx q35[0],q35[1];
h q35[0];
s q35[0];
h q35[0];
s q35[0];
h q35[1];
s q35[1];
cx q35[0],q35[1];
sdg q35[0];
h q35[0];
y q35[0];
s q35[0];
h q35[1];
s q35[1];
z q35[1];
s q35[1];
h q35[1];
s q35[1];
cx q35[0],q35[1];
sdg q35[0];
h q35[0];
sdg q35[1];
h q35[1];
cx q35[0],q35[1];
h q35[0];
s q35[0];
x q35[0];
s q35[0];
h q35[0];
s q35[0];
h q35[1];
s q35[1];
h q35[1];
s q35[1];
cx q35[0],q35[1];
h q35[0];
s q35[0];
sdg q35[1];
h q35[1];
cx q35[0],q35[1];
x q35[0];
s q35[0];
h q35[0];
s q35[0];
x q35[0];
s q35[0];
h q35[0];
s q35[0];
h q35[1];
s q35[1];
z q35[1];
z q35[1];
s q35[1];
cx q35[0],q35[1];
sdg q35[0];
h q35[0];
sdg q35[1];
h q35[1];
cx q35[0],q35[1];
h q35[0];
s q35[0];
h q35[1];
s q35[1];
y q35[1];
h q35[1];
s q35[1];
cx q35[0],q35[1];
sdg q35[0];
h q35[0];
h q35[1];
s q35[1];
cx q35[0],q35[1];
h q35[0];
s q35[0];
x q35[0];
h q35[0];
x q35[1];
h q35[1];
s q35[1];
cx q35[0],q35[1];
x q35[0];
sdg q35[1];
h q35[1];
cx q35[0],q35[1];
sdg q35[0];
h q35[0];
sdg q35[1];
h q35[1];
cx q35[0],q35[1];
x q35[0];
s q35[1];
h q35[1];
s q35[1];
cx q35[0],q35[1];
sdg q35[0];
h q35[0];
z q35[0];
s q35[0];
h q35[1];
s q35[1];
y q35[1];
cx q35[0],q35[1];
h q35[0];
s q35[1];
h q35[1];
s q35[1];
cx q35[0],q35[1];
sdg q35[0];
h q35[0];
x q35[0];
z q35[1];
s q35[1];
cx q35[0],q35[1];
sdg q35[0];
h q35[0];
h q35[1];
s q35[1];
cx q35[0],q35[1];
h q35[0];
z q35[1];
s q35[1];
h q35[1];
s q35[1];
cx q35[0],q35[1];
h q35[0];
s q35[0];
y q35[0];
h q35[0];
s q35[0];
h q35[1];
s q35[1];
z q35[1];
s q35[1];
h q35[1];
s q35[1];
cx q35[0],q35[1];
x q35[0];
x q35[0];
sdg q35[0];
h q35[0];
z q35[1];
h q35[1];
s q35[1];
y q35[1];
s q35[1];
h q35[1];
s q35[1];
cx q35[0],q35[1];
sdg q35[0];
h q35[0];
z q35[0];
sdg q35[0];
h q35[0];
cx q35[0],q35[1];
h q35[0];
s q35[0];
sdg q35[1];
h q35[1];
cx q35[0],q35[1];
z q35[0];
sdg q35[0];
h q35[0];
z q35[1];
cx q35[0],q35[1];
x q35[0];
s q35[0];
h q35[0];
s q35[0];
sdg q35[1];
h q35[1];
cx q35[0],q35[1];
h q35[0];
s q35[0];
x q35[0];
sdg q35[0];
h q35[0];
h q35[1];
s q35[1];
x q35[1];
cx q35[0],q35[1];
sdg q35[0];
h q35[0];
sdg q35[1];
h q35[1];
cx q35[0],q35[1];
h q35[0];
s q35[0];
h q35[1];
s q35[1];
cx q35[0],q35[1];
s q35[0];
h q35[0];
s q35[0];
y q35[1];
h q35[1];
cx q35[0],q35[1];
x q35[0];
h q35[0];
x q35[1];
s q35[1];
cx q35[0],q35[1];
s q35[0];
h q35[0];
s q35[0];
h q35[1];
s q35[1];
x q35[1];
cx q35[0],q35[1];
sdg q35[0];
h q35[0];
z q35[0];
s q35[0];
h q35[0];
s q35[0];
sdg q35[1];
h q35[1];
x q35[1];
s q35[1];
h q35[1];
s q35[1];
cx q35[0],q35[1];
h q35[0];
s q35[0];
z q35[0];
sdg q35[0];
h q35[0];
sdg q35[1];
h q35[1];
y q35[1];
s q35[1];
h q35[1];
s q35[1];
cx q35[0],q35[1];
y q35[0];
s q35[0];
y q35[1];
sdg q35[1];
h q35[1];
cx q35[0],q35[1];
z q35[0];
h q35[0];
s q35[0];
z q35[0];
y q35[1];
sdg q35[1];
h q35[1];
sdg q35[1];
h q35[1];
cx q35[0],q35[1];
sdg q35[0];
h q35[0];
sdg q35[1];
h q35[1];
cx q35[0],q35[1];
s q35[0];
h q35[0];
s q35[0];
x q35[1];
sdg q35[1];
h q35[1];
cx q35[0],q35[1];
y q35[0];
sdg q35[0];
h q35[0];
sdg q35[1];
h q35[1];
y q35[1];
h q35[1];
cx q35[0],q35[1];
h q35[0];
s q35[0];
sdg q35[1];
h q35[1];
cx q35[0],q35[1];
y q35[0];
s q35[0];
s q35[1];
h q35[1];
s q35[1];
cx q35[0],q35[1];
h q35[0];
s q35[0];
h q35[1];
s q35[1];
cx q35[0],q35[1];
x q35[0];
s q35[0];
h q35[0];
s q35[0];
x q35[1];
h q35[1];
cx q35[0],q35[1];
sdg q35[0];
h q35[0];
sdg q35[1];
h q35[1];
cx q35[0],q35[1];
h q35[0];
s q35[0];
h q35[0];
x q35[1];
sdg q35[1];
h q35[1];
cx q35[0],q35[1];
sdg q35[0];
h q35[0];
sdg q35[1];
h q35[1];
cx q35[0],q35[1];
x q35[0];
sdg q35[1];
h q35[1];
cx q35[0],q35[1];
h q35[0];
s q35[0];
h q35[0];
sdg q35[1];
h q35[1];
s q35[1];
h q35[1];
s q35[1];
cx q35[0],q35[1];
sdg q35[0];
h q35[0];
sdg q35[1];
h q35[1];
cx q35[0],q35[1];
h q35[0];
s q35[0];
h q35[1];
s q35[1];
cx q35[0],q35[1];
z q35[0];
sdg q35[0];
h q35[0];
x q35[1];
cx q35[0],q35[1];
sdg q35[0];
h q35[0];
sdg q35[1];
h q35[1];
cx q35[0],q35[1];
h q35[0];
s q35[0];
y q35[0];
s q35[0];
h q35[1];
s q35[1];
z q35[1];
h q35[1];
cx q35[0],q35[1];
sdg q35[0];
h q35[0];
h q35[1];
s q35[1];
cx q35[0],q35[1];
h q35[0];
s q35[0];
x q35[0];
s q35[0];
h q35[0];
s q35[0];
x q35[1];
h q35[1];
cx q35[0],q35[1];
sdg q35[0];
h q35[0];
h q35[1];
s q35[1];
cx q35[0],q35[1];
z q35[0];
s q35[0];
h q35[0];
s q35[0];
h q35[1];
cx q35[0],q35[1];
sdg q35[0];
h q35[0];
sdg q35[1];
h q35[1];
cx q35[0],q35[1];
h q35[0];
s q35[0];
h q35[1];
s q35[1];
cx q35[0],q35[1];
x q35[0];
s q35[0];
h q35[0];
s q35[0];
y q35[1];
s q35[1];
cx q35[0],q35[1];
sdg q35[0];
h q35[0];
h q35[1];
s q35[1];
cx q35[0],q35[1];
h q35[0];
s q35[0];
x q35[0];
s q35[0];
z q35[1];
s q35[1];
h q35[1];
s q35[1];
cx q35[0],q35[1];
sdg q35[0];
h q35[0];
sdg q35[1];
h q35[1];
cx q35[0],q35[1];
h q35[0];
s q35[0];
h q35[1];
s q35[1];
cx q35[0],q35[1];
z q35[0];
s q35[0];
s q35[1];
h q35[1];
s q35[1];
cx q35[0],q35[1];
h q35[0];
s q35[0];
z q35[0];
h q35[0];
z q35[1];
h q35[1];
s q35[1];
cx q35[0],q35[1];
sdg q35[0];
h q35[0];
y q35[0];
h q35[0];
s q35[0];
h q35[1];
cx q35[0],q35[1];
sdg q35[0];
h q35[0];
z q35[0];
h q35[0];
sdg q35[1];
h q35[1];
y q35[1];
sdg q35[1];
h q35[1];
cx q35[0],q35[1];
sdg q35[0];
h q35[0];
sdg q35[1];
h q35[1];
cx q35[0],q35[1];
z q35[0];
sdg q35[0];
h q35[0];
x q35[0];
s q35[0];
h q35[1];
s q35[1];
z q35[1];
s q35[1];
x q35[1];
sdg q35[1];
h q35[1];
cx q35[0],q35[1];
h q35[0];
s q35[0];
z q35[0];
sdg q35[1];
h q35[1];
sdg q35[1];
h q35[1];
cx q35[0],q35[1];
sdg q35[0];
h q35[0];
sdg q35[1];
h q35[1];
cx q35[0],q35[1];
s q35[0];
h q35[0];
s q35[0];
x q35[1];
s q35[1];
h q35[1];
s q35[1];
cx q35[0],q35[1];
sdg q35[0];
h q35[0];
sdg q35[1];
h q35[1];
cx q35[0],q35[1];
h q35[0];
s q35[0];
y q35[0];
s q35[0];
h q35[0];
s q35[0];
h q35[1];
s q35[1];
z q35[1];
sdg q35[1];
h q35[1];
cx q35[0],q35[1];
sdg q35[0];
h q35[0];
h q35[1];
s q35[1];
cx q35[0],q35[1];
y q35[0];
s q35[0];
h q35[0];
s q35[0];
x q35[1];
s q35[1];
cx q35[0],q35[1];
y q35[0];
z q35[1];
s q35[1];
cx q35[0],q35[1];
sdg q35[0];
h q35[0];
sdg q35[1];
h q35[1];
cx q35[0],q35[1];
h q35[0];
s q35[0];
h q35[1];
s q35[1];
cx q35[0],q35[1];
y q35[0];
z q35[1];
s q35[1];
h q35[1];
s q35[1];
cx q35[0],q35[1];
sdg q35[0];
h q35[0];
sdg q35[1];
h q35[1];
cx q35[0],q35[1];
h q35[0];
s q35[0];
x q35[0];
s q35[0];
h q35[1];
s q35[1];
cx q35[0],q35[1];
sdg q35[0];
h q35[0];
y q35[0];
h q35[0];
h q35[1];
s q35[1];
sdg q35[1];
h q35[1];
cx q35[0],q35[1];
h q35[0];
s q35[0];
x q35[0];
h q35[0];
h q35[1];
s q35[1];
z q35[1];
s q35[1];
h q35[1];
s q35[1];
cx q35[0],q35[1];
h q35[0];
s q35[0];
sdg q35[1];
h q35[1];
cx q35[0],q35[1];
y q35[0];
s q35[0];
h q35[1];
s q35[1];
y q35[1];
sdg q35[1];
h q35[1];
cx q35[0],q35[1];
h q35[0];
s q35[0];
z q35[0];
h q35[0];
s q35[0];
sdg q35[1];
h q35[1];
h q35[1];
s q35[1];
cx q35[0],q35[1];
z q35[0];
sdg q35[0];
h q35[0];
h q35[1];
s q35[1];
sdg q35[1];
h q35[1];
cx q35[0],q35[1];
y q35[0];
h q35[0];
s q35[0];
sdg q35[1];
h q35[1];
y q35[1];
cx q35[0],q35[1];
h q35[0];
s q35[0];
h q35[1];
s q35[1];
cx q35[0],q35[1];
x q35[0];
s q35[0];
h q35[0];
s q35[0];
y q35[1];
s q35[1];
h q35[1];
s q35[1];
cx q35[0],q35[1];
h q35[0];
s q35[0];
y q35[0];
s q35[0];
h q35[0];
s q35[0];
sdg q35[1];
h q35[1];
x q35[1];
sdg q35[1];
h q35[1];
cx q35[0],q35[1];
sdg q35[0];
h q35[0];
sdg q35[1];
h q35[1];
cx q35[0],q35[1];
h q35[0];
s q35[0];
x q35[0];
h q35[1];
s q35[1];
y q35[1];
h q35[1];
cx q35[0],q35[1];
sdg q35[0];
h q35[0];
sdg q35[1];
h q35[1];
cx q35[0],q35[1];
h q35[0];
s q35[0];
h q35[1];
s q35[1];
cx q35[0],q35[1];
z q35[0];
h q35[0];
s q35[0];
z q35[1];
sdg q35[1];
h q35[1];
cx q35[0],q35[1];
sdg q35[0];
h q35[0];
x q35[0];
s q35[0];
h q35[0];
s q35[0];
sdg q35[1];
h q35[1];
x q35[1];
sdg q35[1];
h q35[1];
cx q35[0],q35[1];
sdg q35[0];
h q35[0];
sdg q35[1];
h q35[1];
cx q35[0],q35[1];
h q35[0];
s q35[0];
y q35[0];
s q35[0];
y q35[1];
x q35[1];
h q35[1];
cx q35[0],q35[1];
h q35[0];
s q35[0];
sdg q35[1];
h q35[1];
cx q35[0],q35[1];
y q35[0];
h q35[0];
s q35[0];
h q35[1];
s q35[1];
cx q35[0],q35[1];
h q35[0];
s q35[0];
x q35[0];
sdg q35[0];
h q35[0];
z q35[1];
sdg q35[1];
h q35[1];
cx q35[0],q35[1];
z q35[0];
h q35[0];
s q35[0];
h q35[1];
s q35[1];
h q35[1];
s q35[1];
cx q35[0],q35[1];
h q35[0];
s q35[0];
sdg q35[1];
h q35[1];
cx q35[0],q35[1];
x q35[0];
sdg q35[0];
h q35[0];
h q35[1];
s q35[1];
h q35[1];
s q35[1];
cx q35[0],q35[1];
h q35[0];
s q35[0];
h q35[1];
s q35[1];
cx q35[0],q35[1];
x q35[0];
s q35[0];
h q35[0];
s q35[0];
sdg q35[1];
h q35[1];
cx q35[0],q35[1];
h q35[0];
s q35[0];
sdg q35[1];
h q35[1];
cx q35[0],q35[1];
z q35[0];
s q35[0];
h q35[0];
s q35[0];
y q35[1];
h q35[1];
cx q35[0],q35[1];
h q35[0];
s q35[0];
sdg q35[1];
h q35[1];
cx q35[0],q35[1];
x q35[0];
sdg q35[0];
h q35[0];
x q35[1];
s q35[1];
h q35[1];
s q35[1];
cx q35[0],q35[1];
x q35[0];
x q35[1];
sdg q35[1];
h q35[1];
cx q35[0],q35[1];
sdg q35[0];
h q35[0];
z q35[0];
s q35[0];
h q35[1];
s q35[1];
cx q35[0],q35[1];
sdg q35[0];
h q35[0];
h q35[1];
s q35[1];
x q35[1];
cx q35[0],q35[1];
sdg q35[0];
h q35[0];
y q35[0];
sdg q35[0];
h q35[0];
h q35[1];
s q35[1];
h q35[1];
cx q35[0],q35[1];
sdg q35[0];
h q35[0];
sdg q35[1];
h q35[1];
cx q35[0],q35[1];
h q35[0];
s q35[0];
z q35[0];
z q35[1];
h q35[1];
cx q35[0],q35[1];
h q35[0];
s q35[0];
sdg q35[1];
h q35[1];
cx q35[0],q35[1];
x q35[0];
sdg q35[0];
h q35[0];
s q35[1];
cx q35[0],q35[1];
sdg q35[0];
h q35[0];
x q35[0];
s q35[0];
h q35[0];
s q35[0];
h q35[1];
s q35[1];
z q35[1];
h q35[1];
s q35[1];
cx q35[0],q35[1];
y q35[0];
h q35[1];
s q35[1];
x q35[1];
h q35[1];
s q35[1];
cx q35[0],q35[1];
h q35[0];
s q35[0];
s q35[0];
h q35[0];
s q35[0];
y q35[0];
h q35[0];
h q35[1];
s q35[1];
z q35[1];
z q35[1];
h q35[1];
cx q35[0],q35[1];
sdg q35[0];
h q35[0];
s q35[0];
y q35[1];
h q35[1];
cx q35[0],q35[1];
h q35[0];
s q35[0];
x q35[0];
h q35[0];
s q35[0];
h q35[1];
s q35[1];
s q35[1];
cx q35[0],q35[1];
sdg q35[0];
h q35[0];
h q35[0];
s q35[0];
h q35[1];
s q35[1];
y q35[1];
sdg q35[1];
h q35[1];
cx q35[0],q35[1];
sdg q35[0];
h q35[0];
y q35[0];
s q35[0];
h q35[0];
s q35[0];
h q35[1];
s q35[1];
z q35[1];
sdg q35[1];
h q35[1];
cx q35[0],q35[1];
h q35[0];
s q35[0];
sdg q35[1];
h q35[1];
cx q35[0],q35[1];
y q35[0];
h q35[0];
z q35[1];
cx q35[0],q35[1];
sdg q35[0];
h q35[0];
y q35[0];
sdg q35[0];
h q35[0];
h q35[1];
s q35[1];
x q35[1];
s q35[1];
h q35[1];
s q35[1];
cx q35[0],q35[1];
h q35[0];
s q35[0];
z q35[0];
h q35[0];
s q35[0];
sdg q35[1];
h q35[1];
z q35[1];
cx q35[0],q35[1];
h q35[0];
s q35[0];
h q35[0];
s q35[0];
x q35[1];
s q35[1];
h q35[1];
s q35[1];
cx q35[0],q35[1];
sdg q35[0];
h q35[0];
sdg q35[1];
h q35[1];
