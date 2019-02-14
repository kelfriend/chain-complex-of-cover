LoadPackage("HAP");
############################################
#Input:
#-Y: a regular CW-complex
#-G: the fundamental group of Y ("nosimplify")
#-H: a finite index subgroup of G
############################################
#Output:
#-The chain complex of the [G:H]-fold cover of Y, upon which a discrete vector field has been applied
############################################
############################################
ChainComplexOfCover:=function(Y, G, H)
    local trans, gamma, ind, C, basis,
          deform, bij, basedim, DeformCellSgn, nrCriticalCells,
          zero, n, Record, CoverDimension, Boundary;

    if not IsHapRegularCWComplex(Y)
        then
        Error("the provided complex is not a regular CW-complex.\n");
    fi;

    if not IsSubgroup(G,H)
        then
        Error("the given group is not a subgroup of the fundamental group\nof the given CW-complex.\n");
    fi;

    trans:=RightTransversal(G,H);

#C provides the components: basis*, deform, bij and DeformCellSgn*
#*NEEDED TO MODIFY basicRegular.gi TO HAVE THESE

    gamma:=G!.edgeToWord;
    ind:=Index(G,H);
    C:=ChainComplexOfRegularCWComplexWithVectorField(Y);
    basis:=C!.basis;
    deform:=C!.deform;
    bij:=C!.bij;
    basedim:=EvaluateProperty(Y,"dimension");
    DeformCellSgn:=C!.DeformCellSgn;

#Simply returns the number of critical n-cells
############################################
    nrCriticalCells:=function(n);
        return Length(Filtered(Y!.criticalCells, x->x[1]=n));
    end;
############################################

    zero:=[];
    for n in [1..basedim+1]
        do
        zero[n]:=List([1..nrCriticalCells(n-1)],i->0);
    od;

#The list of all previously computed boundaries in the cover
    Record:=List([1..basedim+1],i->[]);

#Gives the number of n-cells in the finite index cover of Y
############################################
    CoverDimension:=function(n)
        return ind*nrCriticalCells(n);
    end;
############################################

#Finds the boundary of the kth n-cell as its image in C_n-1
############################################
    Boundary:=function(n,k)
        local edge, d, b, B, sn, i;

        if n=0
            then
            return [];
        fi;

        if IsBound(Record[n+1][AbsInt(k)])
            then
            return SignInt(k)*Record[n+1][AbsInt(k)];
        fi;

        edge:=gamma(basis[n+1][((k-1) mod nrCriticalCells(n))+1]);
        if k>nrCriticalCells(n)
            then
            d:=Int(Floor(Float(k/nrCriticalCells(n))));
        else
            d:=0;
        fi;

        b:=StructuralCopy(zero[n]);
        B:=Y!.boundaries[n+1][basis[n+1][((k-1) mod nrCriticalCells(n))+1]];
        B:=B{[2..Length(B)]};
        sn:=Y!.orientation[n+1][basis[n+1][((k-1) mod nrCriticalCells(n))+1]];
        B:=List([1..Length(B)],i->sn[i]*B[i]);
        Apply(B,x->DeformCellSgn(n-1,x));
        B:=Concatenation(B);
        Apply(B,i->SignInt(i)*bij[n][AbsInt(i)]);

        for i in B
            do
            b[AbsInt(i)]:=b[AbsInt(i)]+SignInt(i);
        od;

        Record[n+1][k]:=b;
        return Record[n+1][k];
    end;
############################################

    return Objectify(HapChainComplex,
                rec(
                dimension:=CoverDimension,
                boundary:=Boundary,
                properties:=
                [["length", basedim],
                ["connected", true],
                ["type", "chainComplex"],
                ["characteristic",
                0] ]));

end;
############################################
############################################
#Test complexes
2simplices:=
[[1,2,5], [2,5,8], [2,3,8], [3,8,9], [1,3,9], [1,4,9],
[4,5,8], [4,6,8], [6,8,9], [6,7,9], [4,7,9], [4,5,7],
[1,4,6], [1,2,6], [2,6,7], [2,3,7], [3,5,7], [1,3,5]];;
K:=SimplicialComplex(2simplices);;
K:=RegularCWComplex(K);;
Y:=DirectProduct(K,K);
