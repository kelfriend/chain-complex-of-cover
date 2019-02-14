LoadPackage("HAP");
############################################
# Input:
# Y: a regular CW-complex
# G: the fundamental group of Y ("nosimplify")
# H: a finite index subgroup of G
#
# Output:
# The chain complex of the [G:H]-fold cover of Y,
# upon which a discrete vector field has been applied
############################################
############################################
ChainComplexOfCover:=function(Y,G,H)
    local gamma, ind, C, basis, bij,
		  basedim, deform, nrCriticalCells, project, #zero,
		  n, Record, CoverDimension, #RC, RT,
		  MyCanonicalRightCosetElement, deformPair, reindex, canonise, pair2int,
		  SparseToVector, int2pair, VectorToSparse, deindex, PrimaryBoundary,
		  Boundary;

    if not IsHapRegularCWComplex(Y)
        then
        Error("the provided complex is not a regular CW-complex.\n");
    fi;

    if not IsSubgroup(G,H)
        then
        Error("the given group is not a subgroup of the fundamental group\nof the given CW-complex.\n");
    fi;

# C provides the components: deform, bij and DeformCellSgn*
# *NEEDED TO MODIFY basicRegular.gi TO HAVE THIS

    gamma:=G!.edgeToWord;
    ind:=Index(G,H);
    C:=ChainComplexOfRegularCWComplexWithVectorField(Y);
    basis:=C!.basis;
    bij:=C!.bij;
    basedim:=EvaluateProperty(Y,"dimension");
    deform:=C!.DeformCellSgn;

############################################
    nrCriticalCells:=function(n)
        return Length(Filtered(Y!.criticalCells, x->x[1]=n));
    end;
############################################

############################################
	project:=function(n,z)
		return ((z-1) mod nrCriticalCells(n))+1;
	end;
############################################

    #zero:=[];
    #for n in [1..basedim+1]
    #    do
    #    zero[n]:=List([1..ind*nrCriticalCells(n-1)],x->0);
    #od;

# The list of all previously computed boundaries in the cover
    Record:=List([1..basedim+1],i->[]);

# Gives the number of n-cells in the finite index cover of Y
############################################
    CoverDimension:=function(n)
        return ind*nrCriticalCells(n);
    end;
############################################

	#RC:=RightCosets(G,H);
	#RT:=RightTransversal(G,H);
	#RT:=List([1..Length(RT)], x->RT[x]);

############################################
	MyCanonicalRightCosetElement:=function(g)
		local RC, RT, i;
	
		RC:=RightCosets(G,H);
		RT:=RightTransversal(G,H);

		for i in [1..Length(RC)]
			do
			if g in RC[i]
				then
				return RT[i];
			fi;
		od;
	end;
############################################

# The below 8 functions are needed to convert from 'sparse' notation:
# v = [[k,g],[k',g'],[k'',g''],...] (k: a cell, g: a group element)
# to 'vector' notation which will appear as the output:
# v = [0,1,0,0,2,...,2,8,9] (each integer's position uniquely associates it to a pair)
# a reindexing occurs in terms of the critical cells

# takes a sparsely presented pair [k,g] and deforms the cell component
# as per the discrete vector field
############################################
	deformPair:=function(n,p)
		local l, d, j;

		l:=[];
		d:=deform(n,p[1]);
		if not d=[]
			then
			p[1]:=deform(n,p[1]);
			for j in [1..Length(p[1])]
				do
				Add(l, [p[1][j],p[2]]);
			od;
			return l;
		else
			return p;
		fi;
	end;
############################################
 
# reindexes the sparse vector by numbering the critical cells
############################################
	reindex:=function(n,l)
		local i;

		for i in l
			do
			i[1]:=Position(Filtered(Y!.criticalCells, x->x[1]=n), [n,i[1]]);
		od;
	end;
############################################

# makes canonical all group elements in the sparse vector
############################################
	canonise:=function(l)
		local i;

		for i in l
			do
			i[2]:=MyCanonicalRightCosetElement(i[2]);
		od;
	end;
############################################

############################################
	pair2int:=function(p)
		return SignInt(p[1])*(AbsInt(p[1]) + Position(RT, p[2]));
	end;
############################################

############################################
	SparseToVector:=function(n,l)
		local v, i, p;

		Apply(l, x->deformPair(n-1,x));
		reindex(n-1,l);
		canonise(l);

		v:=List([1..CoverDimension(n-1)],x->0);
		for i in l
			do
			p:=pair2int(i);
			v[AbsInt(p)]:=v[AbsInt(p)] + p;
		od;
		return v;
	end;
############################################

############################################
	int2pair:=function(n,z)
		local multiple;

		multiple:=Int((AbsInt(z)-1)/nrCriticalCells(n))+1;
		return [SignInt(z)*project(n,AbsInt(z)),RT[multiple]];
	end;
############################################

############################################
	VectorToSparse:=function(n,l)
		local sparse, i;

		sparse:=[];
		for i in Filtered(l, x-> not x=0)
			do
			Add(sparse, int2pair(n,i));
		od;
		return sparse;
	end;
############################################

# returns the original indexing of the zth critical n-cell
############################################
	deindex:=function(n,z)
		return [n,Filtered(Y!.criticalCells, x->x[1]=n)[z][2]];
	end;
############################################

############################################
	PrimaryBoundary:=function(n,k)
		local pair, b;

		if n=0
			then
			return [];
		fi;

		pair:=int2pair(n,k);
		if n=1
			then
			b:=Y!.boundaries[n+1][deindex(n,AbsInt(pair[1]))[2]];
			return [[SignInt(pair[1])*b[2],RT[1]],[SignInt(pair[1])*b[3],gamma(b[2])]];
		else
			return fail;
		fi;
	end;
############################################

# Finds the boundary of the kth n-cell as its image in C_n-1
############################################
    Boundary:=function(n,k)
        local pb, pair, j; 
			  #b, B, sn, i;

        if n=0
            then
            return [];
        fi;

        if IsBound(Record[n+1][AbsInt(k)])
            then
            return SignInt(k)*Record[n+1][AbsInt(k)];
        fi;

		pb:=PrimaryBoundary(n,k);
		pair:=int2pair(n,k);

		if n=1
			then
			for j in pb
				do
				j[2]:=pair[2]*j[2];
			od;
			
			Record[n+1][AbsInt(k)]:= SparseToVector(n,pb);
        	#b:=StructuralCopy(zero[n]);
        	#B:=Y!.boundaries[n+1][basis[n+1][((k-1) mod nrCriticalCells(n))+1]];
        	#B:=B{[2..Length(B)]};
        	#sn:=Y!.orientation[n+1][basis[n+1][((k-1) mod nrCriticalCells(n))+1]];
        	#B:=List([1..Length(B)],i->sn[i]*B[i]);
        	#Apply(B,x->deform(n-1,x));
        	#B:=Concatenation(B);
        	#Apply(B,i->SignInt(i)*bij[n][AbsInt(i)]);

        	#for i in B
        	#    do
        	#    b[AbsInt(i)]:=b[AbsInt(i)]+SignInt(i);
        	#od;

        	#Record[n+1][k]:=b;

        	#return Record[n+1][k];
	
		else
			return fail;
		fi;
		return Record[n+1][AbsInt(k)];
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
# Test complexes
2simplices:=
[[1,2,5], [2,5,8], [2,3,8], [3,8,9], [1,3,9], [1,4,9],
[4,5,8], [4,6,8], [6,8,9], [6,7,9], [4,7,9], [4,5,7],
[1,4,6], [1,2,6], [2,6,7], [2,3,7], [3,5,7], [1,3,5]];;
K:=SimplicialComplex(2simplices);;
K:=RegularCWComplex(K);;
G:=FundamentalGroupOfRegularCWComplex(K,"n");;
H:=LowIndexSubgroupsFpGroup(G,3)[3];;
GeneratorsOfGroup(H);;
#Y:=DirectProduct(K,K);
