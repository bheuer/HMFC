Attach("Hecke.m");
import "Hecke.m": HeckeOperatorC,HeckeSlopes,HilbertCuspFormCharacter;

F:=QuadraticField(5);
ZF := Integers(F);

for n in [4,8] do
	N:=n*ZF.1*ZF;
	for C in Elements(DirichletGroup(N)) do
		for i1 in [1..4] do
			for i2 in [i1..4] do
				k1:=i1*2+1;
				k2:=i2*2+1;
				
				M:=HilbertCuspFormCharacter(F,N,[k1,k2],C);
				if M`Dimension gt 0 then
				
					PrintFile("job_out4",N);
					PrintFile("job_out4",C);
					PrintFile("job_out4",<k1,k2>);
					PrintFile("job_out4",M`Dimension);
					
					P2 := 2*ZF.1*ZF;F;
					PP2:= Factorisation(P2)[1][1];
					T2:=HeckeOperatorC(M,PP2);
					PrintFile("job_out4",HeckeSlopes(M,PP2));
				end if;
			end for;
		end for;
	end for;
end for;
