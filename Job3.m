Attach("Hecke.m");
import "Hecke.m": HeckeOperatorC,HeckeSlopes,HilbertCuspFormCharacter;

Q<X>:=PolynomialRing(RationalField());

//F<x>:=NumberField(X^4-5*X^2-2*X+1);
F:=QuadraticField(19);
ZF := Integers(F);

for n in [4,8,12] do
	N:=n*ZF.1*ZF;
	for C in Elements(DirichletGroup(N)) do
		for weight in [2*i+1 : i in [1]] do
			M:=HilbertCuspFormCharacter(F,N,[weight,weight],C);
			if M`Dimension gt 0 then
			
				PrintFile("job_out3",N);
				PrintFile("job_out3",C);
				PrintFile("job_out3",<weight,weight>);
				PrintFile("job_out3",M`Dimension);
				
				P2 := 2*ZF.1*ZF;F;
				PP2:= Factorisation(P2)[1][1];
				T2:=HeckeOperatorC(M,PP2);
				PrintFile("job_out3",HeckeSlopes(M,PP2));
			end if;
		end for;
	end for;
end for;
