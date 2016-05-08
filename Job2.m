Attach("Hecke.m");
import "Hecke.m": HeckeOperatorC,HeckeSlopes,HilbertCuspFormCharacter;

F:=QuadraticField(5);
ZF := Integers(F);

for n in [8] do
	N:=n*ZF.1*ZF;
	for C in Elements(DirichletGroup(N)) do
		for weight in [13,15] do
			weight;
			M:=HilbertCuspFormCharacter(F,N,[weight,weight],C);
			if M`Dimension gt 0 then
				PrintFile("job_out2_2",N);
				PrintFile("job_out2_2",C);
				PrintFile("job_out2_2",<weight,weight>);
				PrintFile("job_out2_2",M`Dimension);
				
				P2 := 2*ZF.1*ZF;F;
				PP2:= Factorisation(P2)[1][1];
				T2:=HeckeOperatorC(M,PP2);
				
				//assert T2*T3 eq T3*T2;
				PrintFile("job_out2_2",HeckeSlopes(M,PP2));
			end if;
		end for;
	end for;
end for;
