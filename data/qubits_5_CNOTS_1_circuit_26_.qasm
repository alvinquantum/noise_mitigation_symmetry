OPENQASM 2.0;
include "qelib1.inc";
qreg q106[5];
qreg q107[1];
h q106[3];
cx q107[0],q106[3];
h q106[3];
cx q107[0],q106[2];
cx q107[0],q106[0];
barrier q106[0],q106[1],q106[2],q106[3],q106[4],q107[0];
h q106[0];
s q106[1];
s q106[2];
rz(3.5566106) q106[3];
rz(0.048189056) q106[4];
s q106[4];
s q106[0];
h q106[4];
rz(4.3012379) q106[1];
cx q106[1],q106[0];
s q106[4];
rz(3.4259677) q106[1];
rz(3.867923) q106[4];
barrier q106[0],q106[1],q106[2],q106[3],q106[4],q107[0];
h q106[3];
cx q107[0],q106[3];
h q106[3];
sdg q106[2];
cx q107[0],q106[2];
s q106[2];
h q106[1];
cx q107[0],q106[1];
h q106[1];
h q106[0];
cx q107[0],q106[0];
h q106[0];
