OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
h q[2];
h q[3];
cx q[2],q[0];
cx q[3],q[0];
h q[1];
h q[3];
cx q[3],q[1];
h q[3];
cx q[3],q[2];
h q[3];
rz(1.0*-0.11913647536991812) q[3];
rz(1.0*0.04696493200603326) q[2];
cx q[3],q[2];
rz(1.0*-0.014088325859801982) q[2];
rz(1.0*-0.04051919696006243) q[1];
cx q[3],q[1];
rz(1.0*0.025540610482199102) q[1];
cx q[2],q[1];
rz(1.0*-0.012176400954820124) q[1];
cx q[3],q[1];
rz(1.0*-0.016960316052904708) q[1];
rz(1.0*0.04994656724908739) q[0];
cx q[3],q[0];
rz(1.0*-0.013550325415433915) q[0];
cx q[2],q[0];
rz(1.0*0.025459547948579447) q[0];
cx q[1],q[0];
rz(1.0*0.060765963243375215) q[0];
cx q[3],q[0];
rz(1.0*-0.058502777989906246) q[0];
cx q[2],q[0];
rz(1.0*0.012968433197668605) q[0];
cx q[3],q[0];
rz(1.0*0.04251574032647237) q[0];
cx q[3],q[0];
cx q[1],q[0];
cx q[2],q[1];
cx q[3],q[2];
h q[3];
cx q[3],q[2];
h q[3];
cx q[3],q[1];
h q[3];
h q[1];
cx q[3],q[0];
cx q[2],q[0];
h q[3];
h q[2];
barrier q[0],q[1],q[2],q[3];
h q[0];
h q[2];
cx q[2],q[0];
s q[2];
h q[2];
h q[3];
cx q[3],q[2];
s q[1];
h q[1];
h q[3];
cx q[3],q[1];
s q[3];
h q[3];
rz(1.0*0.007176580991460388) q[3];
rz(1.0*-0.0487976316671623) q[2];
rz(1.0*-0.11553422613600944) q[1];
cx q[3],q[1];
rz(1.0*0.014239700993970934) q[1];
cx q[2],q[1];
rz(1.0*-0.054411439612600976) q[1];
cx q[3],q[1];
rz(1.0*0.013182506426384237) q[1];
rz(1.0*0.010655440441034201) q[0];
cx q[3],q[0];
rz(1.0*-0.08080989761481736) q[0];
cx q[2],q[0];
rz(1.0*0.006551620898334697) q[0];
cx q[3],q[0];
rz(1.0*-0.09900328940640055) q[0];
cx q[1],q[0];
rz(1.0*-0.0955624845963124) q[0];
cx q[3],q[0];
rz(1.0*0.0151112645563178) q[0];
cx q[2],q[0];
rz(1.0*-0.017809415317659633) q[0];
cx q[3],q[0];
rz(1.0*-0.0737525531055947) q[0];
cx q[1],q[0];
cx q[2],q[1];
h q[3];
sdg q[3];
cx q[3],q[1];
h q[3];
h q[1];
sdg q[1];
cx q[3],q[2];
h q[3];
h q[2];
sdg q[2];
cx q[2],q[0];
h q[2];
h q[0];
barrier q[0],q[1],q[2],q[3];
s q[0];
h q[0];
h q[2];
cx q[2],q[0];
h q[2];
s q[1];
h q[1];
h q[3];
cx q[3],q[1];
h q[3];
rz(1.0*0.03062856983384747) q[3];
rz(1.0*0.05901577818895968) q[2];
rz(1.0*-0.0377318703410906) q[1];
cx q[3],q[1];
rz(1.0*0.07053523689139543) q[1];
cx q[2],q[1];
rz(1.0*-0.1108957164449127) q[1];
cx q[3],q[1];
rz(1.0*-0.015478019210044418) q[1];
rz(1.0*-0.03331702630521215) q[0];
cx q[3],q[0];
rz(1.0*-0.002760381534088435) q[0];
cx q[2],q[0];
rz(1.0*0.06692132448957532) q[0];
cx q[3],q[0];
rz(1.0*0.06429688107050259) q[0];
cx q[1],q[0];
rz(1.0*0.04364401596514624) q[0];
cx q[2],q[0];
cx q[3],q[0];
rz(1.0*-0.009715257652999207) q[0];
cx q[3],q[0];
cx q[1],q[0];
cx q[2],q[1];
h q[3];
cx q[3],q[1];
h q[3];
h q[1];
sdg q[1];
h q[2];
cx q[2],q[0];
h q[2];
h q[0];
sdg q[0];
barrier q[0],q[1],q[2],q[3];
h q[1];
h q[2];
cx q[1],q[0];
cx q[2],q[0];
h q[1];
h q[3];
cx q[3],q[1];
h q[3];
cx q[3],q[2];
h q[3];
rz(1.0*-0.041340206405549376) q[3];
rz(1.0*-0.06335401112020661) q[2];
rz(1.0*-0.03352469821449395) q[1];
cx q[3],q[1];
rz(1.0*0.016520423163983333) q[1];
cx q[2],q[1];
rz(1.0*0.0708962003168786) q[1];
rz(1.0*0.008862508246400377) q[0];
cx q[2],q[0];
rz(1.0*-0.01882173642860669) q[0];
cx q[3],q[0];
rz(1.0*-0.047672189605722746) q[0];
cx q[1],q[0];
cx q[3],q[0];
rz(1.0*-0.04056028894547144) q[0];
cx q[2],q[0];
cx q[3],q[0];
rz(1.0*-0.018506486231191492) q[0];
cx q[3],q[0];
cx q[1],q[0];
cx q[2],q[1];
cx q[3],q[1];
h q[3];
cx q[3],q[2];
h q[3];
cx q[3],q[1];
h q[3];
h q[1];
cx q[2],q[0];
cx q[1],q[0];
h q[2];
h q[1];
barrier q[0],q[1],q[2],q[3];
s q[1];
h q[1];
s q[3];
h q[3];
cx q[3],q[1];
s q[2];
h q[2];
s q[3];
h q[3];
cx q[3],q[2];
s q[3];
h q[3];
cx q[3],q[0];
s q[3];
h q[3];
rz(1.0*0.024247692366007434) q[3];
rz(1.0*0.07780226729271573) q[2];
cx q[3],q[2];
rz(1.0*0.059387086859829796) q[2];
rz(1.0*0.057296020835963406) q[1];
cx q[3],q[1];
rz(1.0*0.02518770685820336) q[1];
cx q[2],q[1];
rz(1.0*0.08974416337755287) q[1];
cx q[3],q[0];
rz(1.0*-0.021482116864460638) q[0];
cx q[2],q[0];
rz(1.0*0.036356997764383704) q[0];
cx q[3],q[0];
rz(1.0*0.08831368914017174) q[0];
cx q[1],q[0];
rz(1.0*-0.027917724174571872) q[0];
cx q[3],q[0];
rz(1.0*-0.002548705628870586) q[0];
cx q[2],q[0];
cx q[3],q[0];
rz(1.0*-0.0056813725127048955) q[0];
cx q[1],q[0];
cx q[2],q[1];
cx q[3],q[1];
cx q[3],q[2];
h q[3];
sdg q[3];
cx q[3],q[0];
h q[3];
sdg q[3];
cx q[3],q[2];
h q[3];
sdg q[3];
h q[2];
sdg q[2];
cx q[3],q[1];
h q[3];
sdg q[3];
h q[1];
sdg q[1];
barrier q[0],q[1],q[2],q[3];
s q[0];
h q[0];
h q[1];
h q[2];
h q[3];
cx q[1],q[0];
cx q[2],q[0];
cx q[3],q[0];
h q[1];
h q[2];
h q[3];
rz(1.0*0.037321044683861035) q[0];
cx q[3],q[0];
rz(1.0*-0.020767277613382904) q[0];
cx q[2],q[0];
rz(1.0*0.0787226697415641) q[0];
cx q[3],q[0];
rz(1.0*0.030477151593986705) q[0];
cx q[1],q[0];
rz(1.0*0.06642252150587843) q[0];
cx q[3],q[0];
rz(1.0*-0.041776880309817216) q[0];
cx q[2],q[0];
rz(1.0*0.029648276773298256) q[0];
cx q[3],q[0];
rz(1.0*-0.021458188950601424) q[0];
cx q[1],q[0];
h q[3];
h q[2];
h q[1];
cx q[3],q[0];
cx q[2],q[0];
cx q[1],q[0];
h q[3];
h q[2];
h q[1];
h q[0];
sdg q[0];
barrier q[0],q[1],q[2],q[3];
