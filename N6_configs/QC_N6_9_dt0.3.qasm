OPENQASM 2.0;
include "qelib1.inc";
qreg q[3];
cx q[2],q[1];
h q[2];
rz(0.3*-0.10574811579303006) q[1];
cx q[2],q[0];
rz(0.3*-0.010019939159623972) q[0];
cx q[1],q[0];
rz(0.3*-0.07816253162205503) q[0];
cx q[1],q[0];
cx q[2],q[0];
h q[2];
cx q[2],q[1];
h q[1];
s q[2];
h q[2];
cx q[2],q[0];
s q[2];
h q[2];
rz(0.3*-0.13858711592163908) q[2];
cx q[1],q[0];
rz(0.3*-0.08796557922226136) q[0];
cx q[2],q[0];
rz(0.3*-0.010576831696316874) q[0];
cx q[2],q[0];
cx q[1],q[0];
h q[2];
sdg q[2];
cx q[2],q[0];
h q[2];
sdg q[2];
h q[1];
s q[0];
h q[0];
h q[2];
cx q[2],q[0];
h q[2];
rz(0.3*0.16557277050522973) q[2];
rz(0.3*-0.1892025166162884) q[0];
cx q[2],q[0];
rz(0.3*0.012464299514796745) q[0];
cx q[2],q[0];
h q[2];
cx q[2],q[0];
h q[2];
h q[0];
sdg q[0];
cx q[1],q[0];
h q[1];
cx q[2],q[1];
rz(0.3*-0.09732886577874433) q[1];
rz(0.3*-0.09162232143282573) q[0];
cx q[1],q[0];
rz(0.3*-0.026834130053148914) q[0];
cx q[1],q[0];
cx q[2],q[1];
h q[1];
cx q[1],q[0];
s q[1];
h q[1];
h q[2];
cx q[2],q[0];
h q[2];
cx q[2],q[1];
rz(0.3*-0.17686269592902545) q[1];
cx q[2],q[0];
rz(0.3*0.0031319523146368655) q[0];
cx q[1],q[0];
rz(0.3*-0.08901933254017533) q[0];
cx q[1],q[0];
cx q[2],q[0];
cx q[2],q[1];
h q[2];
cx q[2],q[0];
h q[2];
h q[1];
sdg q[1];
