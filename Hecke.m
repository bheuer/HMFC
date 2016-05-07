
Attach("/usr/local_machine/magma/package/Geometry/ModFrmHil/definite.m");
import "/usr/local_machine/magma/package/Geometry/ModFrmHil/definite.m": InvariantSpace;

Attach("/usr/local_machine/magma/package/Geometry/ModFrmHil/precompute.m");
import "/usr/local_machine/magma/package/Geometry/ModFrmHil/precompute.m": get_tps;


function findScaling(PL,N,v1,v2)
	//return c such that c*v1 = v2
	
	P1Rep:=PL`P1Rep;
	_,_,f1:=P1Rep(v1,false,true);
	_,_,f2:=P1Rep(v2,false,true);
	
	//now we have f1*v1 = f2*v2
	return f1*Modinv(f2,N);
end function;	

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
    //was chi:=[(Modinv(x[1][1],N)*r[1][1]) mod N : r in res ];

    chi:=[findScaling(ProjLine,N,x,r) : r in res ];
    

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
                    bool, u0, a := P1Rep(u,true,true); //a*u=u0
                    if bool then
						// at this point we have a*u=u0;
						
                        elt_data:=lookup[u0];
                        n:=Index(CFDm, elt_data[1]);
                        if n ne 0 then //is this a contributing orbit? If not, don't need to compute
                            
                            x_n := FundamentalDomain[elt_data[1]]; //==FundamentalDomain[CFDm[n]] by def of index n
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
                           
                            T:= WR(fac)*(C(c));
                            
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

function HilbertCuspFormCharacter(F,N,weight,C)
	
	M:=HilbertCuspForms(F,N,weight);
	_:=Dimension(M);//computes everything
	
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

