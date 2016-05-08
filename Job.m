Attach("Hecke.m");
import "Hecke.m": HeckeOperatorC,HeckeSlopes,HilbertCuspFormCharacter;

Q<X>:=PolynomialRing(RationalField());

//F<x>:=NumberField(X^4-5*X^2-2*X+1);
F<x>:=QuadraticField(3);
ZF<g> := Integers(F);

for n in [4,8,12] do
	N:=n *ZF;
	for C in Elements(DirichletGroup(N)) do
		for weight in [2*i+1 : i in [1..3]] do
			M:=HilbertCuspFormCharacter(F,N,[weight,weight],C);
			if M`Dimension gt 0 then
			
				//PrintFile("job_out",C);
				//PrintFile("job_out",<weight,weight>);
				//PrintFile("job_out",M`Dimension);
				
				P2 := 2*ZF.1*ZF;F;
				PP2:= Factorisation(P2)[1][1];
				T2:=HeckeOperatorC(M,PP2);
				P3 := 3*ZF.1*ZF;
				PP3:= Factorisation(P3)[1][1];
				T3:=HeckeOperatorC(M,PP3);
				
				assert T2*T3 eq T3*T2;
				HeckeSlopes(M,PP2);
				//PrintFile("job_out",HeckeSlopes(M,PP2));
			end if;
		end for;
	end for;
end for;
