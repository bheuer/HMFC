
Attach("/usr/local_machine/magma/package/Geometry/ModFrmHil/definite.m");
import "/usr/local_machine/magma/package/Geometry/ModFrmHil/definite.m": InvariantSpace;

Attach("/usr/local_machine/magma/package/Geometry/ModFrmHil/precompute.m");
import "/usr/local_machine/magma/package/Geometry/ModFrmHil/precompute.m": get_tps;


function twisted_invariant_space(M,C,weight_field,alpha)

    O:=M`QuaternionOrder;
    B:=Algebra(O);

    F := BaseField(B);
    ZF:= Integers(F);
    N := M`Level;

    HMDFs := M`ModFrmHilDirFacts; 
    hmdf:=HMDFs[1];
    ProjLine:=hmdf`PLD; 
    FundamentalDomain:=ProjLine`FD;
    SplittingMap:=ProjLine`splitting_map;
    Stabs:=ProjLine`Stabs;
    P1Rep:=ProjLine`P1Rep;
    
    WR := hmdf`weight_rep;
	weight_dim := hmdf`weight_dimension;
	
	weight_field2:=weight_field;
	//weight_field2:=Compositum(weight_field,Codomain(C));//TODO:Does that actually work?
	M2K:=MatrixRing(weight_field2, weight_dim);
	WR_ := map<B-> M2K | q:-> WR(q)>;


	S := [B| u@Umap : u in Generators(U)] where U, Umap := UnitGroup(BaseField(B));
	S cat:= [s[1] : s in Stabs[alpha] | s[1] notin {B|1,-1}];

	Gamma:= [SplittingMap(s) : s in S];
	x:=FundamentalDomain[alpha];

	res := [g*x : g in Gamma];
	
	if x[1][1] ne 0 and Gcd(Norm(N),Norm(x[1][1]*ZF)) eq 1 then
		chi:=[(Modinv(x[1][1],N)*r[1][1]) mod N : r in res ];
	else
		chi:=[(Modinv(x[2][1],N)*r[2][1]) mod N : r in res ];
	end if;
	WR_chi_1:=map<B -> M2K|q :-> WR_(q)*(C(chi[Position(S,q)])^(-1))>;
	
	L := InvariantSpace(Stabs[alpha],WR_chi_1,weight_dim,weight_field2);
	return L;
end function;

function basis_matrix(M,C)

    HMDFs := M`ModFrmHilDirFacts; 
    hmdf:=HMDFs[1];
    ProjLine:=hmdf`PLD; 
    stabs:=ProjLine`Stabs;
    weight_field:= hmdf`weight_base_field;
    weight_field2:=Compositum(weight_field,Codomain(C));

	l:=1;
	repeat
	  L:=twisted_invariant_space(M,C,weight_field2,l);
	  if Rank(L) eq 0 then l:=l+1; end if;
	until (l gt #stabs) or (Rank(L) ne 0);
	
	if l gt #stabs then
	  return [],0,0;
	end if;
	contrib_orbs:=[l];
	for m0:=l+1 to #stabs do
		N:=twisted_invariant_space(M,C,weight_field2,m0);
		if Rank(N) ne 0 then
			Append(~contrib_orbs, m0);
			nb_rows:=Nrows(L)+Nrows(N);
			nb_cols:=Ncols(L)+Ncols(N);
			Q:=RMatrixSpace(weight_field2, nb_rows, nb_cols)!0;
			InsertBlock(~Q, L, 1, 1);
			InsertBlock(~Q, N, Nrows(L)+1, Ncols(L)+1);
			L:=Q;
		end if;
	end for;
	contrib_orbs := {@ x : x in contrib_orbs @};
	N:=Transpose(Solution(Transpose(L), ScalarMatrix(weight_field2, Rank(L), 1)));
	
	return contrib_orbs, L, N;
end function;

function compute(M,PP,C,CFD)
	O:=M`QuaternionOrder;
    B:=Algebra(O);
    F:=BaseField(B);
    ZF:=Integers(F);
    N:=M`Level;
    assert Domain(C) eq F;

    HMDFs := M`ModFrmHilDirFacts;
    hmdf:=HMDFs[1];
    ProjLine:=hmdf`PLD;
    FundamentalDomain:=ProjLine`FD;
    SplittingMap:=ProjLine`splitting_map;
    Stabs:=ProjLine`Stabs;
    P1Rep:=ProjLine`P1Rep;
    lookup:=ProjLine`Lookuptable;
    max_order_units:=hmdf`max_order_units;
    //assume Cl O = 1 for now
    Tps:=get_tps(M,PP*ZF.1);
	
    //m=1
    CFDm:=CFD;//hmdf`CFD //contributing fundamental domain
    CFDl:=CFD;//hmdf`CFD;

    WR := hmdf`weight_rep;
    wd := hmdf`weight_dimension;

    weight_field:= hmdf`weight_base_field;
    weight_field2:=Compositum(weight_field,Codomain(C));
    
    Tp := Matrix(weight_field2, wd*#CFDl, wd*#CFDm, []);

    for l :=1 to #HMDFs do
        for m:=1 to #HMDFs do
			ts:=Tps[<m,l>];
			
			for ll:=1 to #ts do
				for mm:=1 to #CFDm do //basically sum over fundamental domain, but only contributing orbits matter
					mat:=SplittingMap(ts[ll]);
				
					x_m :=FundamentalDomain[mm];
					u := mat*x_m;
					bool, u0 := P1Rep(u,true,false);
					if bool then
						elt_data:=lookup[u0];
						n:=Index(CFDm, elt_data[1]);
						if n ne 0 then //is this a contributing orbit? If not, don't need to compute
							x_n := FundamentalDomain[elt_data[1]]; //==FundamentalDomain[CFDm[n]] by def of index n
							o:=max_order_units[elt_data[2]];
							
							//At this point we have
							//x_mm pi^{-1} ^O_0*  = t^{-1}o x_n  ^O_0*
							//where o is the maximal order unit that                 
							//translates the rep in the fundamental domain           
							//x_n to the element (t x_mm pi^{-1} ^O_0*) in           
							//P^1(ZF/N). So we will later mult out o^{-1}*t

							fac := o^(-1)*ts[ll];
							//finally, we have to adjust pi as follows:
							// we know that it is in ^O_0*, but we want it
							// to be in ^O_1*. So we just need to multiply
							// by an appropriate element of Z(A) and then
							// mult out by chi( this element ). In order to          
							// do so, note that we know
							//           (d  *)
							//      pi = (0  *)
							// for some invertible element d at pp. So we            
							// just need to find d. To this end, reorder as
							//
							//    o^{-1} t x_mm = x_n ^O_1* u*pi
							//
							// Then we can just compute the upper left entry
							// of o^{-1} t x_mm (this does not dep on                
							// choice of x_mm) and the upper left entry of
							// x_n (does not dep on choice of x_n) and thus
							// we get
							//
							//   a(o^{-1} t x_mm) = a(x_n) *a(u*pi)
							//
							// This gives a:=a(u*pi) mod PP. We then                 
							// just need to multiply through by a to get
							// a suitable representative pi' = a^{-1}u*pi
							// which is in O_1 and then we get
							//
							// x_m pi'^{-1} ^O_1*= t^{-1}o x_n ^O_1*a^{-1}
							//
							// Therefore, we compute the desired image to be
							//
							// f(x_m pi^{-1}) = f(x_n)^{o^{-1}t} chi(a)^{-1}
							
							a1 := (SplittingMap(o^(-1))*SplittingMap(ts[ll])*x_m)[1][1];
							a2 := x_n[1][1];//this is an element of the FD so should hopefully be invertible at least after reiteration
							a:= a1*a2^(-1);
							
							M:= WR(fac)*(C(a))^(-1); //although this is C(a^(-1)),
													 //it is easier to just invert C(a)
													 //instead of finding a modular
													 //inverse of a mod PP
							// note that WR(gamma) is -^{gamma}, so this is the right element to act with
							// we now have the matrix for which f(x_m pi^{-1}) = f(x_n) * M

							//So we just need to add everything up:
							X := ExtractBlock(Tp, (n-1)*wd+1, (mm-1)*wd+1, wd, wd);
							InsertBlock(~Tp, X+M, (n-1)*wd+1, (mm-1)*wd+1);
						end if;
					end if;
				end for;
			end for;
        end for;
    end for;

return Tp;

end function;


//This is the Hecke operator one should actually use. Insert invariants (basis matrix) for bm
function HeckeOperatorC(M,PP,C,bm,bminv,CFD)
    T:=compute(M,PP,C,CFD);
    bmT := bm * ChangeRing(T, BaseRing(bm));
    TM := bmT * bminv;
    return TM;
end function;

function HeckeEigenvalues(T,PP)
    f:=CharacteristicPolynomial(Matrix(T));
    return Slopes(NewtonPolygon(f,PP));
end function;




