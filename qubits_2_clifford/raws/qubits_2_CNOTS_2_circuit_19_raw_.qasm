OPENQASM 2.0;
include "qelib1.inc";
qreg q18[2];
s q18[0];
sdg q18[1];
h q18[1];
cx q18[0],q18[1];
h q18[0];
s q18[0];
sdg q18[1];
h q18[1];
cx q18[0],q18[1];
