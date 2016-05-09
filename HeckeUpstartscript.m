Attach("Hecke.m");
import "Hecke.m": HeckeOperatorC,HeckeSlopes,HilbertCuspFormCharacter;

F<x>:=QuadraticField(5);
ZF<g> := Integers(F);
PP2:= Factorisation(2*ZF.1*ZF)[1][1];
PP3:= Factorisation(3*ZF.1*ZF)[1][1];
N:= 4*ZF;

C:=[i : i in Elements(DirichletGroup(N))][2];

for k in [2*i+1: i in [1..7]] do
	M:=HilbertCuspFormCharacter(F,N,[k,k],C);
		
	T2:=HeckeOperatorC(M,PP2);
	HeckeSlopes(M,PP2);
end for;

//T3:=HeckeOperatorC(M,PP3);
//assert T2*T3 eq T3*T2;
