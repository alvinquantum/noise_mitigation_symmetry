OPENQASM 2.0;
include "qelib1.inc";
qreg q41[2];
s q41[0];
h q41[0];
s q41[0];
cx q41[0],q41[1];
sdg q41[0];
h q41[0];
h q41[0];
s q41[0];
y q41[1];
cx q41[0],q41[1];
h q41[0];
s q41[0];
sdg q41[0];
h q41[0];
sdg q41[1];
h q41[1];
cx q41[0],q41[1];
sdg q41[0];
h q41[0];
s q41[0];
sdg q41[1];
h q41[1];
cx q41[0],q41[1];
sdg q41[0];
h q41[0];
x q41[0];
sdg q41[1];
h q41[1];
x q41[1];
s q41[1];
h q41[1];
s q41[1];
cx q41[0],q41[1];
h q41[0];
s q41[0];
sdg q41[1];
h q41[1];
cx q41[0],q41[1];
z q41[0];
h q41[1];
s q41[1];
x q41[1];
h q41[1];
s q41[1];
cx q41[0],q41[1];
sdg q41[0];
h q41[0];
x q41[0];
s q41[0];
h q41[0];
s q41[0];
sdg q41[1];
h q41[1];
z q41[1];
cx q41[0],q41[1];
h q41[0];
s q41[0];
sdg q41[1];
h q41[1];
