Attach("Hecke.m");
import "Hecke.m": HeckeOperatorC,HeckeSlopes,HilbertCuspFormCharacter;

F<x>:=QuadraticField(7);
ZF<g> := Integers(F);
PP3:= Factorisation(3*ZF.1*ZF)[1][1];
PP2:= Factorisation(2*ZF.1*ZF)[1][1];
N:= 4*ZF;

for C in [i : i in Elements(DirichletGroup(N))] do
	M:=HilbertCuspFormCharacter(F,N,[5,5],C);
	if M`Dimension eq 0 then continue; end if;
	M;
	T2:=HeckeOperatorC(M,PP2);
	HeckeSlopes(M,PP2);
end for;
