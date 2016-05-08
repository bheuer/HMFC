Attach("Hecke.m");
import "Hecke.m": HeckeOperatorC,HeckeSlopes,HilbertCuspFormCharacter;

F<x>:=QuadraticField(5);
ZF<g> := Integers(F);
P2 := 2*ZF.1*ZF;
PP2:= Factorisation(P2)[1][1];
N:= 8*ZF;

for C in Elements(DirichletGroup(N)) do
	M:=HilbertCuspFormCharacter(F,N,[9,9],C);
	if M`Dimension eq 0 then
		continue;
	end if;
	
	M;

	T2:=HeckeOperatorC(M,PP2);
	
	HeckeSlopes(M,PP2);
end for;

