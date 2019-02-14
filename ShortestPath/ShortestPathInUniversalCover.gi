LoadPackage("HAP");
#LoadPackage("kbmag");
LoadPackage("nq");
#######################################################
#######################################################
ShortestPath:=function(Y,G,vS,g)
	local elts, nqelts, visited, distance, skeleton, strictskeleton, phi, theta, Neighbours, Update, c, H, w, x, q, j, k, altdistance, p, pp;

    if not IsHapRegularCWComplex(Y) then
        Print("The given complex is not a regular CW-complex.\n");
        return fail;
    fi;

	if not g in G then
		Print("The given group element does not belong to the fundamental group of the given CW-complex.\n");
		return fail;
	fi;

#Creates a list of some elements of the fundamental group of Y
    elts:=[1..Y!.nrCells(0)];; Apply(elts, x->[One(G),g]);

#Assigns a weight of 1 to all edges if weighting is not already present
    if not IsBound(Y!.weights) then
		Print("The given complex is unweighted. Uniform weighting has been applied.\n");
        Y!.weights:=[1..Y!.nrCells(1)]; Apply(Y!.weights, x -> 1);
    fi;

#Initialisation step: all distances from source are set to infinity (bar vS) and all cells are marked as unvisited.
    visited:=[1..Y!.nrCells(0)]; Apply(visited, x->[0,0]);
    distance:=[1..Y!.nrCells(0)]; Apply(distance, x->[infinity,infinity]); distance[vS][1]:=0;

#Extracts the 1-skeleton from Y
	skeleton:=ShallowCopy(Y!.boundaries[2]);; Apply(skeleton, x->[x[2],x[3]]); Apply(skeleton, x->Set(x));
	strictskeleton:=ShallowCopy(Y!.boundaries[2]);; Apply(strictskeleton, x->[x[2],x[3]]);

#For non-abelian groups, creates an epimorphism from G to a nilpotent-quotient subgroup
	phi:=IdentityMapping(G);
	if not IsAbelian(G) then
		phi:=NqEpimorphismNilpotentQuotient(G,10);
	fi;

#Creates a second list of group elements under the nq epimorphism to ease up on calculation
	nqelts:=List(elts,x->[Image(phi,x[1]),Image(phi,x[2])]);

#Knuth-Bendix?
	#theta:=KBMAGRewritingSystem(G);
	#SetInfoLevel(InfoRWS, 1);
	#KnuthBendix(theta);

#Finds the neighbours of a given vertex via the edgeToWord function
#######################################################
    Neighbours:=function(v,h)
		local neighbours, bone, i;
        neighbours:=[];
        for bone in skeleton do
            if v in bone then
                Add(neighbours, Filtered(bone, x->not x=v)[1]);
            fi;
        od;
        for i in neighbours do
            if [v,i] in strictskeleton then
                neighbours[Position(neighbours,i)]:=[i,h*G!.edgeToWord(Position(strictskeleton, [v,i]))];
            else
                neighbours[Position(neighbours,i)]:=[i,h*(G!.edgeToWord(Position(strictskeleton, [i,v])))^-1];
            fi;
        od;
        return neighbours;
    end;
#######################################################

#Ensures that each relevant group element is stored in elts (& nqelts) and that entries in 'distance' and 'visited' stay defined
#######################################################
    Update:=function(vv,hh)
        local N, n, x, l, k, u;
        N:=Neighbours(vv,hh);
        for n in N do
             x:=Image(phi,n[2]);
             if not x in nqelts[n[1]] then
                Add(elts[n[1]],n[2]);
                Add(nqelts[n[1]],x);
                l:=Length(elts[n[1]]);
#should a group element not be in elts, it is added and an extra Y!.nrCells(0) entries are
#added to both 'distance' and 'visited' to account for the added layer in the universal cover
                distance[n[1]][l]:=infinity;
                visited[n[1]][l]:=0;
            fi;
        od;
        return N;
    end;
#######################################################

#Main body; Dijkstra's algorithm for shortest path in a weighted graph
	c:=[vS,One(G)];
    	p:=1;
	if g=One(G) then
		return 0;	#this corresponds to a trivial loop
	else
    	while visited[vS][2]=0 do
        	q:=Update(c[1],c[2]);
        	for j in q do
#the algorithm cycles through the neighbours of the selected vertex, c,
#and redefines the distance to each neighbour of c if the weight applied
#to the edge connecting c to its neighbour plus the distance from vS to c
#be less than the existing entry in 'distance' for that neighbour
            	pp:=Position(nqelts[j[1]],Image(phi,j[2]));
        		altdistance:=
        		distance[c[1]][p]
        		+Y!.weights[Position(skeleton, Set([c[1],j[1]]))];
        		if altdistance < distance[j[1]][pp] then
        			distance[j[1]][pp]:=altdistance;
        		fi; 
        	od;
#the previous choice for c is now declared 'visited' and the next choice is
#made by finding the (first) unvisited vertex whose distance from vS is minimal.
			visited[c[1]][p]:=1;			
			H:=[];
			for w in [1..Y!.nrCells(0)] do
				for x in [1..Length(distance[w])] do
					if visited[w][x]=0 then
						Add(H, [distance[w][x],[w,x]]);
					fi;
				od;
			od;
			H:=MinimumList(H);
			p:=H[2][2];
			c:=[H[2][1], elts[H[2][1]][p]];
    	od;
#'distance' was given two entries for each vertex at the initialisation step; one for the
#base space, and one for the layer defined by g
#therefore the output will always lie in the same predetermined place
    	return distance[vS][2];
	fi;
end;
#######################################################
#######################################################
Y:=[[1,1,1],[1,0,1],[1,1,1]];
Y:=PureCubicalComplex(Y);
Y:=RegularCWComplex(Y);

2simplices:=
[[1,2,5], [2,5,8], [2,3,8], [3,8,9], [1,3,9], [1,4,9],
[4,5,8], [4,6,8], [6,8,9], [6,7,9], [4,7,9], [4,5,7],
[1,4,6], [1,2,6], [2,6,7], [2,3,7], [3,5,7], [1,3,5]];;
K:=SimplicialComplex(2simplices);
YY:=RegularCWComplex(K);

G:=FundamentalGroupOfRegularCWComplex(Y, "nosimplify");
GG:=FundamentalGroupOfRegularCWComplex(YY, "nosimplify");
