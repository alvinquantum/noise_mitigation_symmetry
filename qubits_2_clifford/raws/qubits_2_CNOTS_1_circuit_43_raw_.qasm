OPENQASM 2.0;
include "qelib1.inc";
qreg q42[2];
y q42[0];
sdg q42[0];
h q42[0];
z q42[1];
s q42[1];
h q42[1];
s q42[1];
cx q42[0],q42[1];
h q42[0];
s q42[0];
h q42[1];
s q42[1];
