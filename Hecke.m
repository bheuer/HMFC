
Attach("/usr/local_machine/magma/package/Geometry/ModFrmHil/definite.m");
import "/usr/local_machine/magma/package/Geometry/ModFrmHil/definite.m": InvariantSpace;

Attach("/usr/local_machine/magma/package/Geometry/ModFrmHil/precompute.m");
import "/usr/local_machine/magma/package/Geometry/ModFrmHil/precompute.m": get_tps;


function twisted_invariant_space(M,alpha)

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
    
    C:=M`DirichletCharacter;
    weight_field := M`weight_base_field;
    WR := hmdf`weight_rep;
    weight_dim := hmdf`weight_dimension;
    
    M2K:=MatrixRing(weight_field, weight_dim);
    WR_ := map<B-> M2K | q:-> WR(q)>;


    S := [B| u@Umap : u in Generators(U)] where U, Umap := UnitGroup(BaseField(B));
    S cat:= [s[1] : s in Stabs[alpha] | s[1] notin {B|1,-1}];

    Gamma:= [SplittingMap(s) : s in S];
    x:=FundamentalDomain[alpha];

    res := [g*x : g in Gamma];
    
    //in fact, after our update, x[1][1] is always invertible
    chi:=[(Modinv(x[1][1],N)*r[1][1]) mod N : r in res ];
    //if x[1][1] ne 0 and (ZF!1 in x[1][1]*ZF.1*ZF+N) then
    //    
    //else
    //   chi:=[(Modinv(x[2][1],N)*r[2][1]) mod N : r in res ];
    //end if;//one of these must work because P1(ZF/N) has coprime entries
    
    WR_chi_1:=map<B -> M2K|q :-> WR_(q)*(C(chi[Position(S,q)])^(-1))>;
    
    L := InvariantSpace(Stabs[alpha],WR_chi_1,weight_dim,weight_field);
    return L;
end function;

function basis_matrix(M,hmdf)
	
    ProjLine:=hmdf`PLD; 
    stabs:=ProjLine`Stabs;
    weight_field:= hmdf`weight_base_field;
    
    contrib_orbs:=[];
    L:=Matrix(weight_field, 0, 0,[]);
    
    for m0:=1 to #stabs do
        N:=twisted_invariant_space(M,m0);
        if Rank(N) ne 0 then //This orbit contributes non-trivial functions
            Append(~contrib_orbs, m0);
            nb_rows:=Nrows(L)+Nrows(N);
            nb_cols:=Ncols(L)+Ncols(N);
            Q:=RMatrixSpace(weight_field, nb_rows, nb_cols)!0;
            InsertBlock(~Q, L, 1, 1);
            InsertBlock(~Q, N, Nrows(L)+1, Ncols(L)+1);
            L:=Q;
        end if;
    end for;
    contrib_orbs := {@ x : x in contrib_orbs @};
    N:=Transpose(Solution(Transpose(L), ScalarMatrix(weight_field, Rank(L), 1)));
    
    hmdf`CFD := contrib_orbs;
    hmdf`basis_matrix := L;
    hmdf`basis_matrix_inv := N;
    return hmdf;
    
end function;

function computeHeckeMatrix(M,PP)
    O:=M`QuaternionOrder;
    B:=Algebra(O);
    F:=BaseField(B);
    ZF:=Integers(F);
    N:=M`Level;
    C:=M`DirichletCharacter;

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
    CFDm:=hmdf`CFD; //contributing fundamental domain
    CFDl:=hmdf`CFD;

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
                            
                            
							"u";u;
							"u0";u0;
                            
                            x_n := FundamentalDomain[elt_data[1]]; //==FundamentalDomain[CFDm[n]] by def of index n
                            o:=max_order_units[elt_data[2]];
                            
                            lhs:= SplittingMap(ts[ll])*x_m;
                            rhs:=SplittingMap(o)*x_n;
                            
                            "lhs";lhs;"rhs";rhs;
                            
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
                            
                            SplittingMap(o^(-1))*SplittingMap(ts[ll])*x_m;
							
                            a1 := (SplittingMap(o^(-1))*SplittingMap(ts[ll])*x_m)[1][1];
							
                            a2 := x_n[1][1];//this is an element of the FD so should hopefully be invertible at least after reiteration
                            a:= a1*a2^(-1);
                            
                            
							if not ZF!1 in (a*ZF.1*ZF+N) then
								"hazard";
								continue;
							end if;
								
                            
                            T:= WR(fac)*(C(a))^(-1); //although this is C(a^(-1)),
                                                     //it is easier to just invert C(a)
                                                     //instead of finding a modular
                                                     //inverse of a mod PP
                            // note that WR(gamma) is -^{gamma}, so this is the right element to act with
                            // we now have the matrix for which f(x_m pi^{-1}) = f(x_n) * M

                            //So we just need to add everything up:
                            X := ExtractBlock(Tp, (n-1)*wd+1, (mm-1)*wd+1, wd, wd);
                            InsertBlock(~Tp, X+T, (n-1)*wd+1, (mm-1)*wd+1);
                        end if;
                    end if;
                end for;
            end for;
        end for;
    end for;

return Tp;

end function;


//This is the Hecke operator one should actually use. Insert invariants (basis matrix) for bm
function HeckeOperatorC(M,PP)
	
    HMDFs := M`ModFrmHilDirFacts; 
    hmdf:=HMDFs[1];
    bm:=hmdf`basis_matrix;
    bminv:=hmdf`basis_matrix_inv;
    

    T:=computeHeckeMatrix(M,PP);
    bmT := bm * ChangeRing(T, BaseRing(bm));
    TM := bmT * bminv;
    return TM;
end function;

function HeckeEigenvalues(T,PP)
    f:=CharacteristicPolynomial(Matrix(T));
    return Slopes(NewtonPolygon(f,PP));
end function;


procedure adjustFundamentalDomain(M)
	N:=M`Level;
	O:=M`QuaternionOrder;
	B:=Algebra(O);
	F:=BaseField(B);
	ZF:=Integers(F);
	
	HMDF:=M`ModFrmHilDirFacts;
	for i:=1 to #HMDF do
		hmdf := HMDF[i];
		units:=hmdf`max_order_units;
		UU := Universe(units);
		
		ProjLine:=hmdf`PLD;
		
		sm:=ProjLine`splitting_map;
		FD:=ProjLine`FD;
		P1Rep:=ProjLine`P1Rep;
		for xindex:=1 to #FD do
			x:=FD[xindex];
			if  (ZF!1) notin (x[1][1]*ZF+N) then //first coeff of rep is not invertible
				for uindex :=1 to #units do
					u:=units[uindex];
					translate := (sm(u)*x);
					if ZF!1 in translate[1][1]*ZF+N then //if translate invertible
						//good, now we need to update a few thing.
						
						//just translate the stabilisers by inverse
						uinv := u^(-1);

						ProjLine`Stabs[xindex] := [ [* UU!(gamma[1]*uinv), 1*] :gamma in ProjLine`Stabs[xindex]];//nobody knows why, but elements of Stab[i] must first be unmasked from some "*". Hence the [1].
						
						//StabOrders doesn't change of course
						//Lookuptable doesn't either, because we have chosen gamma*x (and not gamma*x*d for some d in (ZF/N)*).
						
						//Most importantly, because that's the point of all of this:
						_,ProjLine`FD[xindex]:=P1Rep(sm(u)*x,false,false);
						//hopefully that's enough
						break;
					end if;
				if uindex ge #units then
					//I really hope you never get here, but I can't prove it right now.

					assert false;
				end if;
				end for;
				
			end if;
		end for;
		hmdf`PLD:=ProjLine;
		M`ModFrmHilDirFacts[i]:=hmdf;//THIS LINE IS NEEDED! #TalkingAboutMemory
	end for;
end procedure;

function HilbertCuspFormCharacter(F,N,weight,C)
	
	M:=HilbertCuspForms(F,N,weight);
	_:=Dimension(M);//computes everything

	//need to tweak projective line until all reps in the
	//fundamental domain have mod invertible first entry.
	adjustFundamentalDomain(M);
	
	is_character_trivial:=(Order(C) eq 1);
	M`DirichletCharacter:=C;
	if not is_character_trivial then//this is to make sure the old Hecke operators still work in the case of trivial character
		M`weight_base_field := Compositum(M`weight_base_field,Codomain(C));
	end if;
	
	HMDF := M`ModFrmHilDirFacts;
	
	dim:=0;
	for i:=1 to #HMDF do
		HMDF[i]`weight_base_field := M`weight_base_field;
	
		hmdf :=basis_matrix(M,HMDF[i]);
		HMDF[i]:=hmdf;
		dim+:=Nrows(hmdf`basis_matrix);
	end for;
	
	if is_character_trivial then
		assert M`Dimension eq dim; //sanity check
	else
		M`Dimension := dim;
	end if;
	
	M`ModFrmHilDirFacts := HMDF; //probably unneccesary, but who nows how magma deal with memory
	
	return M;

end function;

