LoadPackage("HAP");
#############################Input##################################
########### Y: a regular CW-complex ################################
########### G: the fundamental group of Y ("nosimplify") ###########
########### H: a finite index subgroup of G ########################
#############################Output#################################
########### C: the chain complex of the [G:H]-fold cover of Y, #####
############## with a discrete vector field. #######################
####################################################################

ChainComplexOfCover:=function(Y,G,H)

	local
		  gamma, ind, c, basedim, deform, nrCriticalCells,
		  firstOfPair, SparseRecord, Record, CoverDimension,
		  RC, RT, MyCanonicalRightCosetElement, deformPair,
		  reindex, canonise, pair2int, SparseToVector, int2pair,
		  VectorToSparse, deindex, gRemove, rearrange, 
		  SparseBoundary, Boundary;

	if not IsHapRegularCWComplex(Y)
		then
		Error(
		"the provided complex is not a regular CW-complex.\n"
		);
	fi;
	if not IsSubgroup(G,H)
		then
		Error(
		"the given group is not a subgroup of the fundamental group"
		);
	fi;

# basicRegular.gi was altered to yield DeformCellSgn

	gamma:=G!.edgeToWord;
	ind:=Index(G,H);
	c:=ChainComplexOfRegularCWComplexWithVectorField(Y);
	basedim:=EvaluateProperty(Y,"dimension");
	deform:=c!.DeformCellSgn;


# outputs the number of critical n-cells in C_n
####################################################################
	nrCriticalCells:=function(n)
		return Length(Filtered(Y!.criticalCells, x->x[1]=n));
	end;
####################################################################


# takes a reindexed vector and ouputs which critical n-cell it
# is referring to in that indexing
####################################################################
	firstOfPair:=function(n,z)
		return ((z-1) mod nrCriticalCells(n))+1;
	end;
####################################################################


# the dictionaries of all previously computed boundaries
	SparseRecord:=NewDictionary([], true);
	Record:=NewDictionary([], true);


# gives the number of n-cells in the finite index cover of Y
####################################################################
	CoverDimension:=function(n)
		return ind*nrCriticalCells(n);
	end;
####################################################################


	RC:=RightCosets(G,H);
	RT:=RightTransversal(G,H);


# inputs an arbitrary group element and returns the canonical right
# coset element correspoding to it
####################################################################
	MyCanonicalRightCosetElement:=function(g)
		local i;

		for i in [1..Length(RC)]
			do
			if g in RC[i]
				then
				return RT[i];
			fi;
		od;
	end;
####################################################################


####################################################################
############# The below 9 functions are needed to ##################
############## convert v from 'sparse' notation: ###################
####################################################################
################### v = [[k,g],[k',g'],[k'',g''],...] ##############
################### k: the kth n-cell, #############################
################### g: some group element ##########################
####################################################################
############# to 'vector' notation: ################################
####################################################################
################### v = [0,1,0,-2,...,2,8,9] #######################
################### each integer's position in v uniquely ##########
################### associates it to a sparse pair #################
################### note: a reindexing occurs in terms #############
################### of the critical cells ##########################
####################################################################


# takes a sparsely presented pair [k,g] and deforms the cell
# components, outputting a list of all deformed cells all with the
# same group element
####################################################################
	deformPair:=function(n,p)
		local l, d, j;

		l:=[];
		d:=deform(n,p[1]);
		if not d=[]
			then
			for j in [1..Length(d)]
				do
				Add(l, [d[j],p[2]]);
			od;
			return l;
		else
			return [];
		fi;
	end;
####################################################################


# reindexes the sparse vector in terms of the critical cells
####################################################################
	reindex:=function(n,l)
		local i;

		for i in l
			do
			i[1]:=SignInt(i[1])*Position(
				Filtered(Y!.criticalCells,x->x[1]=n),
                [n,AbsInt(i[1])]
				);
		od;
	end;
####################################################################


# makes canonical all group elements in the sparse vector
####################################################################
	canonise:=function(l)
		local i;

		for i in l
			do
			i[2]:=MyCanonicalRightCosetElement(i[2]);
		od;
	end;
####################################################################


# takes a sparse pair and uniquely associates it to an integer z
# thus, the zth position in 'vector' notation corresponds to this
# pair
####################################################################
	pair2int:=function(n,p)
        local pos;

        pos:=Position(RT, p[2]);
		return SignInt(p[1])*(AbsInt(p[1])+
        (pos-1)*nrCriticalCells(n));
	end;
####################################################################


# uses the above 4 functions to finally convert to vector notation
####################################################################
	SparseToVector:=function(n,l)
		local v, i, p;

		Apply(l, x->deformPair(n-1,x));
		l:=Concatenation(l);
		canonise(l);
		reindex(n-1,l);
		v:=List([1..CoverDimension(n-1)],x->0);
		for i in l
			do
			p:=pair2int(n-1,i);
			v[AbsInt(p)]:=v[AbsInt(p)] + SignInt(p);
			od;
		return v;
	end;
####################################################################


# inverse function of pair2int
####################################################################
	int2pair:=function(n,z)
		local multiple;

		multiple:=Int((AbsInt(z)-1)/nrCriticalCells(n))+1;
		return [SignInt(z)*firstOfPair(n,AbsInt(z)),
				RT[multiple]];
	end;
####################################################################


# returns a sparse vector in its canonical, reindexed form
####################################################################
	VectorToSparse:=function(n,l)
		local sparse, i;

		sparse:=[];
		for i in Filtered(l, x-> not x=0)
			do
			Add(sparse, int2pair(n,i));
		od;
		return sparse;
	end;
####################################################################


# returns the original indexing of the zth critical n-cell
####################################################################
	deindex:=function(n,z)
		return [n,SignInt(z)*
				Filtered(Y!.criticalCells,x->x[1]=n)[AbsInt(z)][2]];
	end;
####################################################################


# takes a sparse list and removes its group elements
####################################################################
	gRemove:=function(l)
		local r, cell;

		r:=[];
		for cell in l
			do
			Add(r, cell[1]);
		od;
		return r;
	end;
####################################################################


# rearranges the sparse boundary to match the order given in
# Y!.boundaries so that orientations are consistent
####################################################################
	rearrange:=function(l,b)
		local eltless, r, i;

		eltless:=gRemove(b);
		Apply(eltless, x->AbsInt(x));
		r:=[];
		for i in l
			do
			r[Position(eltless, AbsInt(i[1]))]:=i;
		od;
		return r;
	end;
####################################################################


# inputs a sparse pair as well as its dimension, n, and outputs
# a sparsely presented boundary list of n-1 cells in the cover
####################################################################
	SparseBoundary:=function(n,s)

        local
			  b, o, b1, g, b2, b3, jigsaw,
			  0jigsaw, delta1, 0delta1, tick, delta2,
			  neighbour, 0delta2, inter,
			  h, newG, i, sb;

		if n=0
			then
			return [];
		fi;

		if not LookupDictionary(SparseRecord,[n,s])=fail
			then
			return LookupDictionary(SparseRecord,[n,s]);
		fi;

		b:=Y!.boundaries[n+1][AbsInt(s[1])];
		b:=b{[2..Length(b)]};
		o:=Y!.orientation[n+1][AbsInt(s[1])];

		if n=1
			then
			b1:=[o[1]*SignInt(s[1])*b[1],s[2]];
			g:=gamma(s[1]);
			b2:=[o[2]*SignInt(s[1])*b[2],s[2]*g];
            AddDictionary(SparseRecord,[n,s],[b1,b2]);
		else
			b3:=ShallowCopy(b);
			Apply(b3, x->[SignInt(s[1])*x,s[2]]);
			jigsaw:=[];
			Add(jigsaw, b3[1]);

			delta1:=SparseBoundary(n-1,b3[1]);
			0delta1:=gRemove(delta1);
			Apply(0delta1, x->AbsInt(x));

			tick:=1;
			while not Length(jigsaw)=Length(b3)
				do
				0jigsaw:=gRemove(jigsaw);
				Apply(0jigsaw, x->AbsInt(x));
				for neighbour in Filtered(b3{[tick+1..Length(b3)]},
				x->not x[1] in 0jigsaw)
					do
					delta2:=SparseBoundary(n-1,neighbour);
					0delta2:=gRemove(delta2);
					Apply(0delta2, x->AbsInt(x));

					inter:=Intersection(0delta1,0delta2);
					if not inter=[]
						then
						h:=
						delta2[Position(0delta2, inter[1])][2];
						newG:=
						delta1[Position(0delta1, inter[1])][2];
						Add(jigsaw,
						[neighbour[1],h*newG*neighbour[2]]
						);
					fi;
					if Length(jigsaw)>tick
						then
						break;
					fi;
				od;
				tick:=tick+1;
				if IsBound(jigsaw[tick])
					then
					delta1:=SparseBoundary(n-1,jigsaw[tick]);
					0delta1:=gRemove(delta1);
					Apply(0delta1, x->AbsInt(x));
				fi;
			od;
			for i in [1..Length(jigsaw)]
				do
				jigsaw[i][1]:=jigsaw[i][1]*o[i];
			od;
			sb:=rearrange(jigsaw,b3);
			AddDictionary(SparseRecord,[n,s],sb);
		fi;
		return LookupDictionary(SparseRecord,[n,s]);
	end;
####################################################################


# Finds the boundary of the kth n-cell as its image in C_n-1
####################################################################
	Boundary:=function(n,k)
		local pair, prebound;

		if not LookupDictionary(Record,[n,k])=fail
			then
			return LookupDictionary(Record,[n,k]);
		fi;

		pair:=int2pair(n,k);
        pair[1]:=deindex(n,pair[1])[2];
        prebound:=StructuralCopy(SparseBoundary(n,pair));

        AddDictionary(Record,[n,k],SparseToVector(n,prebound));
		return LookupDictionary(Record, [n,k]);
	end;
####################################################################


	return Objectify(HapChainComplex,
		rec(
			dimension:=CoverDimension,
			deform:=deform,
			sparseBoundary:=SparseBoundary,
			boundary:=Boundary,
			properties:=
				[["length", basedim],
				["connected", true],
				["type", "chainComplex"],
				["characteristic",0]]
			));

end;
####################################################################
####################################################################
# Test complex
K:=PureCubicalKnot(3,1);;
K:=PureComplexComplement(K);;
Y:=RegularCWComplex(K);;
G:=FundamentalGroupOfRegularCWComplex(Y, "no");;
H:=LowIndexSubgroupsFpGroup(G,5);
H:=H[5];