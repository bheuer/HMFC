delete HeckeOperatorC;
delete HeckeSlopes;
delete HilbertCuspFormCharacter;

Attach("Hecke.m");
import "Hecke.m": HeckeOperatorC,HeckeSlopes,HilbertCuspFormCharacter;

F<x>:=QuadraticField(3);
ZF<g> := Integers(F);
N:= 8*ZF;

C:=[c : c in Elements(DirichletGroup(N)) | Order(c) le 2][5];
M:=HilbertCuspFormCharacter(F,N,[5,5],C);

HMDFs := M`ModFrmHilDirFacts;
hmdf:=HMDFs[1];
ProjLine:=hmdf`PLD;
FD:=ProjLine`FD;

P2 := 2*ZF.1*ZF;
PP2:= Factorisation(P2)[1][1];
T2:=HeckeOperatorC(M,PP2);
P3 := 3*ZF.1*ZF;
PP3:= Factorisation(P3)[1][1];
T3:=HeckeOperatorC(M,PP3);
assert T2*T3 eq T3*T2;
