OPENQASM 2.0;
include "qelib1.inc";
qreg q28[2];
y q28[0];
s q28[0];
y q28[1];
h q28[1];
s q28[1];
cx q28[0],q28[1];
h q28[0];
s q28[0];
z q28[0];
h q28[0];
s q28[0];
z q28[1];
h q28[1];
s q28[1];
cx q28[0],q28[1];
x q28[0];
sdg q28[0];
h q28[0];
x q28[0];
s q28[0];
h q28[0];
s q28[0];
x q28[1];
s q28[1];
h q28[1];
s q28[1];
y q28[1];
sdg q28[1];
h q28[1];
sdg q28[1];
h q28[1];
cx q28[0],q28[1];
h q28[0];
s q28[0];
s q28[0];
h q28[0];
s q28[0];
h q28[1];
s q28[1];
z q28[1];
sdg q28[1];
h q28[1];
cx q28[0],q28[1];
sdg q28[0];
h q28[0];
sdg q28[1];
h q28[1];
cx q28[0],q28[1];
z q28[0];
sdg q28[0];
h q28[0];
h q28[1];
s q28[1];
z q28[1];
h q28[1];
s q28[1];
cx q28[0],q28[1];
z q28[0];
s q28[0];
h q28[0];
s q28[0];
sdg q28[1];
h q28[1];
s q28[1];
cx q28[0],q28[1];
h q28[0];
s q28[0];
x q28[1];
s q28[1];
h q28[1];
s q28[1];
cx q28[0],q28[1];
sdg q28[1];
h q28[1];
