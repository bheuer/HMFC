
PP2:= 2*ZF.1*ZF;
PP3:= 3*ZF.1*ZF;
PP7:= 7*ZF.1*ZF;

O:=M`QuaternionOrder;
B:=Algebra(O);
F:=BaseField(B);
ZF:=Integers(F);
N:=M`Level;

HMDFs := M`ModFrmHilDirFacts;
hmdf:=HMDFs[1];
ProjLine:=hmdf`PLD;
FD:=ProjLine`FD;
SplittingMap:=ProjLine`splitting_map;
Stabs:=ProjLine`Stabs;
P1Rep:=ProjLine`P1Rep;
lookup:=ProjLine`Lookuptable;
max_order_units:=hmdf`max_order_units;
weight_dim := hmdf`weight_dimension;
weight_field:= hmdf`weight_base_field;
weight_field2:=Compositum(weight_field,Codomain(C));//TODO:Does that actually work?


T2:=HeckeOperatorC(M,PP2);
T3:=HeckeOperatorC(M,PP3);
T7:=HeckeOperatorC(M,PP7);

if Order(C) eq 1 then
	assert T2 eq HeckeOperator(M,PP2);
	assert T3 eq HeckeOperator(M,PP3);
end if;
assert T2*T3 eq T3*T2;

f:=CharacteristicPolynomial(T2);
NP:=NewtonPolygon(f,PP2);
Slopes(NP);
