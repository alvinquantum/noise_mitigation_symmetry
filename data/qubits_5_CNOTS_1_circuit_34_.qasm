OPENQASM 2.0;
include "qelib1.inc";
qreg q138[5];
qreg q139[1];
h q138[4];
cx q139[0],q138[4];
h q138[4];
cx q139[0],q138[3];
cx q139[0],q138[2];
barrier q138[0],q138[1],q138[2],q138[3],q138[4],q139[0];
s q138[0];
rz(5.4801816) q138[1];
h q138[1];
s q138[4];
rz(1.5626891) q138[0];
h q138[0];
rz(1.752116) q138[1];
cx q138[2],q138[1];
rz(1.8538157) q138[0];
swap q138[4],q138[0];
rz(6.0808313) q138[0];
barrier q138[0],q138[1],q138[2],q138[3],q138[4],q139[0];
cx q139[0],q138[3];
cx q139[0],q138[2];
cx q139[0],q138[1];
h q138[0];
cx q139[0],q138[0];
h q138[0];
