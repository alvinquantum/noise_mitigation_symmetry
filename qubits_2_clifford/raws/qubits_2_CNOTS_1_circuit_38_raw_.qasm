OPENQASM 2.0;
include "qelib1.inc";
qreg q37[2];
z q37[0];
s q37[0];
x q37[1];
sdg q37[1];
h q37[1];
cx q37[0],q37[1];
sdg q37[0];
h q37[0];
sdg q37[1];
h q37[1];
