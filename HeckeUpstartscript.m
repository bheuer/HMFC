delete compute;
delete basis_matrix;
delete HeckeOperatorC;
delete HeckeEigenvalues;
delete twisted_invariant_space;

Attach("Hecke.m");
import "Hecke.m": twisted_invariant_space,compute,HeckeOperatorC,basis_matrix,HeckeEigenvalues;

F<x>:=QuadraticField(5);
ZF := Integers(F);
N:= 8*ZF;

M:=HilbertCuspForms(F,N,[3,3]);
Dimension(M);//computes everything

C:=[c : c in Elements(DirichletGroup(N)) | Order(c) le 2][3];

PP2:= 2*ZF.1*ZF;
PP3:= 3*ZF.1*ZF;
PP7:= 7*ZF.1*ZF;

CFD,L,L_inv:=basis_matrix(M,C);


O:=M`QuaternionOrder;
B:=Algebra(O);
F:=BaseField(B);
ZF:=Integers(F);
N:=M`Level;

HMDFs := M`ModFrmHilDirFacts;
hmdf:=HMDFs[1];
ProjLine:=hmdf`PLD;
FundamentalDomain:=ProjLine`FD;
SplittingMap:=ProjLine`splitting_map;
Stabs:=ProjLine`Stabs;
P1Rep:=ProjLine`P1Rep;
lookup:=ProjLine`Lookuptable;
max_order_units:=hmdf`max_order_units;
weight_dim := hmdf`weight_dimension;
weight_field:= hmdf`weight_base_field;
weight_field2:=Compositum(weight_field,Codomain(C));//TODO:Does that actually work?
	

T2:=HeckeOperatorC(M,PP2,C,L,L_inv,CFD);
T3:=HeckeOperatorC(M,PP3,C,L,L_inv,CFD);
//T7:=HeckeOperatorC(M,PP7,C,L,L_inv,CFD);

assert T2*T3 eq T3*T2;

f:=CharacteristicPolynomial(T2);
NP:=NewtonPolygon(f,PP2);
Slopes(NP);
