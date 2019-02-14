LoadPackage("HAP");
############################################
#Input:
#-Y: a regular CW-complex
#-G: the fundamental group of Y
#-H: a finite index subgroup of G
############################################
#Output:
#-The chain complex of the [G:H]-fold cover of Y, upon which a discrete vector field has been applied
############################################
############################################
ChainComplexOfCover:=function(Y, G, H)
    local ind, C, deform, Record, CoverDimension, Boundary;

    if not IsHapRegularCWComplex(Y)
        then
        Error("the provided complex is not a regular CW-complex.\n");
    fi;

    if not IsSubgroup(G,H)
        then
        Error("the given group is not a subgroup of the fundamental group of the given CW-complex.\n");
    fi;

#C provides a deform function which acts on the discrete vector field of the cover 
    ind:=Index(G,H);
    C:=ChainComplexOfRegularCWComplexWithVectorField(Y);
    deform:=C!.deform;

#The list of all previously computed boundaries in the cover
    Record:=List([1..ind],
            x->[List([1..EvaluateProperty(Y, "dimension")+1], 
            y->[])]);

#Gives the number of n-cells in the finite index cover of Y
############################################
    CoverDimension:=function(n)
        return ind*Length(Filtered(Y!.criticalCells, x->x[1]=n));
    end;
############################################

#Finds the boundary of the jth k-cell as its image in C_k-1 
#(determined by the gth transversal element)
############################################
    Boundary:=function(g,k,j)
        local n, b, l;
        
        if k=0
            then
            return [];
        fi;
            
        n:=j/AbsInt(j);
        j:=AbsInt(j);

        if IsBound(Record[g][k+1][j]) or IsBound(Record[g][k+1][n*j])
            then
            return Apply(Record[g][k+1][j], i->i*n);
        fi;

        if [k,j] in Y!.criticalCells
            then
            b:=[];
            for l in [2..Length(Y!.boundaries[k+1][j])]
                do
                Add(b, Y!.boundaries[k+1][j][l]);
            od;
            Apply(b, x->deform(k-1,x));
            b:=Product(b,Y!.orientation[k+1][j]);
            b:=Concatenation(b);
            Record[g][k+1][j]:=b;
        else
            Error("no such cell exists in the cover.\n");
        fi;

        return Record[g][k+1][j];
    end;
############################################

    return Objectify(HapChainComplex,
                rec(
                dimension:=CoverDimension,
                boundary:=Boundary,
                properties:=
                [["length", Dimension(Y)],
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
