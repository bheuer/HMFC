
Attach("/usr/local_machine/magma/package/Geometry/ModFrmHil/definite.m");
import "/usr/local_machine/magma/package/Geometry/ModFrmHil/definite.m": InvariantSpace;

Attach("/usr/local_machine/magma/package/Geometry/ModFrmHil/precompute.m");
import "/usr/local_machine/magma/package/Geometry/ModFrmHil/precompute.m": get_tps;


function findScaling(PL,N,v1,v2)
	//given two equivalent representatives v1, v2 of
	//the projective line P1(ZF/N),
	//return c such that c*v1 = v2
	
	P1Rep:=PL`P1Rep;
	_,_,f1:=P1Rep(v1,false,true);
	_,_,f2:=P1Rep(v2,false,true);
	
	//now we have f1*v1 = f2*v2
	return f1*Modinv(f2,N);
end function;	

function twisted_invariant_space(M,hmdf,m)
	//m is the fundamental domain index

    O:=M`QuaternionOrder;
    B:=Algebra(O);

    F := BaseField(B);
    ZF:= Integers(F);
    N := M`Level;

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
    S cat:= [s[1] : s in Stabs[m] | s[1] notin {B|1,-1}];

    Gamma:= [SplittingMap(s) : s in S];
    x:=FundamentalDomain[m];

    res := [g*x : g in Gamma];
    
    //in fact, after our update, x[1][1] is always invertible
    //was chi:=[(Modinv(x[1][1],N)*r[1][1]) mod N : r in res ];

    chi:=[findScaling(ProjLine,N,x,r) : r in res ];

    WR_chi_1:=map<B -> M2K|q :-> WR_(q)*(C(chi[Position(S,q)])^(-1))>;
    
    L := InvariantSpace(Stabs[m],WR_chi_1,weight_dim,weight_field);
    return L;
end function;

function basis_matrix(M,hmdf)
	
    ProjLine:=hmdf`PLD; 
    stabs:=ProjLine`Stabs;
    weight_field:= hmdf`weight_base_field;
    
    contrib_orbs:=[];
    L:=Matrix(weight_field, 0, 0,[]);
    
    for m0:=1 to #stabs do
        N:=twisted_invariant_space(M,hmdf,m0);
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


function big_basis_matrix(M)
	//write all separate matrices into one large block matrix
	HMDF := M`ModFrmHilDirFacts;
	nrows := &+ [Nrows(HMDF[m]`basis_matrix): m in [1..#HMDF]];
	ncols := &+ [Ncols(HMDF[m]`basis_matrix): m in [1..#HMDF]];
	B := Matrix(BaseRing(HMDF[1]`basis_matrix), nrows, ncols, []);
	row := 1;
	col := 1;
	for hmdf in HMDF do
		if not IsEmpty(hmdf`CFD) then
			InsertBlock(~B, hmdf`basis_matrix, row, col);
			row +:= Nrows(hmdf`basis_matrix);
			col +:= Ncols(hmdf`basis_matrix);
		end if;
	end for;
	Binv := Transpose(Solution(Transpose(B), IdentityMatrix(BaseRing(B), Nrows(B))));
	
	M`basis_matrix_big := B;
	M`basis_matrix_big_inv := Binv;
	return M;
end function;


function computeHeckeMatrix(M,PP)
	//compute big Hecke matrix

	cached, Tp := IsDefined(M`HeckeBig,PP);
	if cached then
		return Tp;
	end if;
	
    O:=M`QuaternionOrder;
    B:=Algebra(O);
    F:=BaseField(B);
    ZF:=Integers(F);
    N:=M`Level;
    C:=M`DirichletCharacter;
	SplittingMap := M`splitting_map;

	HMDFs := M`ModFrmHilDirFacts;
    
    weight_field:= HMDFs[1]`weight_base_field;
    Tps:=get_tps(M,PP);//I am sorry
    
    HMDFs := M`ModFrmHilDirFacts;
    wd:=M`weight_dimension;
    
    Tp := MatrixRing(weight_field,Ncols(M`basis_matrix_big)) ! 0;//seems to be a nice way to describe the zero matrix: just coerce zero into ring
	
	inds := [l : l in [1..#HMDFs] | #(HMDFs[l]`CFD) ne 0];

	row := 0;
	col := 0;
    for m in inds do
		hmdfm:=HMDFs[m];
		PLm:=hmdfm`PLD;		//Projective LineT2
		FDm:=PLm`FD;		//Fundamental Domain
		CFDm:=hmdfm`CFD;	//fundamental domain reps of contributing orbits
		
        for l in inds do
			defined,ts:= IsDefined(Tps,<m,l>);
            if defined then //don't try to just continue to get rid if the if! We need to increase row in any case
            
				hmdfl	:=HMDFs[l];
				PLl		:=hmdfl`PLD;	//Projective Line
				FDl		:=PLl`FD;		//Fundamental Domain
				lookup	:=PLl`Lookuptable;
				P1Rep	:=PLl`P1Rep;
				CFDl:=hmdfl`CFD;		//fundamental domain reps of contributing orbits
				max_order_units := hmdfl`max_order_units;

				WR := hmdfl`weight_rep;
			
				Tpml:= Matrix(weight_field, wd*#CFDl, wd*#CFDm, []);
				
				for ll:=1 to #ts do
					mat:=SplittingMap(ts[ll]);
					for mm:=1 to #CFDm do //basically sum over fundamental domain, but only contributing orbits matter
						x_m :=FDm[CFDm[mm]];
						u := mat*x_m;
						bool, u0, a := P1Rep(u,true,true); //a*u=u0
						if bool then
							// at this point we have a*u=u0;
							
							elt_data:=lookup[u0];
							n:=Index(CFDl, elt_data[1]);
							if n ne 0 then //is this a contributing orbit? If not, don't need to compute
								
								x_n := FDl[elt_data[1]]; //==FundamentalDomain[CFDl[n]] by def of index n
								o:=max_order_units[elt_data[2]];
								
								_,_,b:=P1Rep(SplittingMap(o)*x_n,false,true);
								//at this point we have 
								//b*o*x_n * ^O_1= u0 = a*t*x_m ^O_1.
								//We therefore set:
								c:=a*Modinv(b,N);
								
								//such that after comparing pi-valuations we get 
								
								//t*x_m *c*pi^(-1) ^O_1* = o*x_n ^O_1*
								
								//where pi is some element with pi-norm 1 st
								//the upper left coefficient is 1, ie it is 
								//in ^O_1 (without units). Moreover,
								//o is the maximal order unit that                 
								//translates the rep in the fundamental domain           
								//x_n to the element (t x_mm pi^{-1} ^O_0*) in           
								//P^1(ZF/N). So we will later mult out o^{-1}*t

								fac := o^(-1)*ts[ll];
								
								// So when we just mult out c via the char C.
								// We thus compute the desired image to be:
								//
								// f(x_m pi^{-1}) = f(x_n)^{o^{-1}t} chi(c)^{-1}
							   
								//TODO: check whether I have to invert chi here.
								
								T:= WR(fac)*C(c);
								
								// note that WR(gamma) is -^{gamma}, so this is the right element to act with
								// we now have the matrix for which f(x_m pi^{-1}) = f(x_n) * M
								
								//So we just need to add everything up:
								X := ExtractBlock(Tpml, (n-1)*wd+1, (mm-1)*wd+1, wd, wd);
								InsertBlock(~Tpml, X+T, (n-1)*wd+1, (mm-1)*wd+1);
							end if;
						end if;
					end for;
				end for;
				InsertBlock(~Tp, Tpml, row+1,col+1);
			end if;
            row +:= (#(HMDFs[l]`CFD)) * wd;
        end for;
        col +:= #(HMDFs[m]`CFD)* wd;
		row := 0;
    end for;
    M`HeckeBig[PP]:=Tp;//cache
return Tp;

end function;


//This is the Hecke operator one should actually use. Insert invariants (basis matrix) for bm
function HeckeOperatorC(M,PP)
	
    HMDFs := M`ModFrmHilDirFacts;
    bm:=M`basis_matrix_big;
    bminv:=M`basis_matrix_big_inv;
    
    T:=computeHeckeMatrix(M,PP);
    bmT := bm * ChangeRing(T, BaseRing(bm));
    TM := bmT * bminv;
    return TM;
end function;

function HeckeSlopes(M,PP)
	T_PP:=HeckeOperatorC(M,PP);//is cached, so hopefully not too not inefficient
    
    //coefficients lie in weight_base_field, not nec. in F
    WF:=M`weight_base_field;
    OWF:=Integers(WF);
    
    PP_OWF:=0*OWF; //==PP*OWF
    for g in Generators(PP) do
		PP_OWF+:=g*OWF;
	end for;
	fact:=Factorisation(PP_OWF)[1]; //doesn't matter which one for valuation because extension is cyclotomic => Galois
	
	PPP:=fact[1];
	e:=fact[2];
	
	f:=CharacteristicPolynomial(Matrix(T_PP));
	NP:=NewtonPolygon(f,PPP);
	S:=[s/e : s in Slopes(NP)];
	V:=AllVertices(NP);	
	
	swm := [];
	lastslope:=Infinity();
	
	for i:=1 to #V-1 do
		v1:=V[i];
		v2:=V[i+1];
		mult:=v2[1]-v1[1];
		slope:=(v1[2]-v2[2])/mult;
		if #swm gt 0 and slope eq lastslope then
			swm[#swm]:=<slope, mult+swm[#swm][2]>;
		else 
			Append(~swm,<slope/e,mult>);
			lastslope:=slope;
		end if;
	end for;
    return swm;
    
end function;

function HilbertCuspFormCharacter(F,N,weight,C)
	
	M:=HilbertCuspForms(F,N,weight);
	_:=Dimension(M);//does all precomputation of MFHDs
	
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
	
	M`ModFrmHilDirFacts := HMDF; //probably unneccesary, but who nows how magma deal with memory //update: yep, necessary
	M:=big_basis_matrix(M);
	
	return M;

end function;

