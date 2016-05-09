Attach("Hecke.m");
import "Hecke.m": HeckeOperatorC,HeckeSlopes,HilbertCuspFormCharacter;

F<x>:=QuadraticField(5);
ZF<g> := Integers(F);
P2 := 2*ZF.1*ZF;
PP2:= Factorisation(P2)[1][1];

ideal11 := 11*ZF;
fact:= Factorisation(ideal11);

data:=fact[1];
PP11:=data[1];
e:=data[2];

N:=8*PP11;
C:=[i : i in Elements(DirichletGroup(8*ZF))| Order(i) eq 2][2];

assert C(1)   eq 1;
assert C(-1)  eq 1;
assert C(g)   eq 1;
assert C(g+4) eq -1;

for C in [i : i in Elements(DirichletGroup(8*ZF))|Order(i) eq 2 ] do
	M:=HilbertCuspFormCharacter(F,N,[4,4],C);
	if M`Dimension gt 0 then
		M;
		HeckeSlopes(M,PP2);	
	end if;
end for;

