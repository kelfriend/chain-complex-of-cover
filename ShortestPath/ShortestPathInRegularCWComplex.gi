LoadPackage("HAP");

##########################################################
##########################################################
ShortestPath:=function(Y,vS,vT)
    local h, unvisited, distance, previous, skeleton, strictskeleton, Neighbours, c, k, altdistance, U, path, properpath, l, m;

    if not IsHapRegularCWComplex(Y) then
        Print("The given complex is not a regular CW-complex.\n");
        return fail;
    fi;

    if not EvaluateProperty(Y, "dimension") < infinity then
        Print("The given complex is not finite.\n");
        return fail;
    fi;

    if not IsBound(Y!.weights) then
        Print("The given complex is unweighted. Uniform weighting has been applied.\n");
        Y!.weights:=[1..Y!.nrCells(1)]; Apply(Y!.weights, x -> 1);
    fi;
    
#Initialisation step: all distances from source are set to infinity (bar vS) and all cells are marked as unvisited
    unvisited:=[1..Y!.nrCells(0)]; Apply(unvisited, x -> 0);
    distance:=[1..Y!.nrCells(0)]; Apply(distance, x -> infinity); distance[vS]:=0;

#This will track the previously visited vertex in the optimal path from the source to that vertex
    previous:=[];

#Gives the 1-skeleton of the CW complex
    skeleton:=ShallowCopy(Y!.boundaries[2]);; Apply(skeleton, x->Set([x[2],x[3]]));
    strictskeleton:=ShallowCopy(Y!.boundaries[2]);; Apply(strictskeleton, x->[x[2],x[3]]);

##########################################################
    Neighbours:=function(v)
        local coboundary, neighbours, i, j;
        coboundary:=[];
        neighbours:=[];
        for i in [2..Length(Y!.coboundaries[1][v])] do
            Add(coboundary, Y!.coboundaries[1][v][i]);            
        od;
        for j in coboundary do;
            Add(neighbours, Filtered(skeleton[j], x -> not x=v)[1]);
        od;
        return neighbours;
    end;
##########################################################

#Obtains shortest distance from vS to vT as well as a list of previously visited vertices in the optimal path
    c:=vS;
    while not c = vT do
        for k in Neighbours(c) do
            altdistance:= distance[c] + Y!.weights[Position(skeleton, Set([c,k]))];
            if altdistance < distance[k] then
                distance[k]:= altdistance;
                previous[k]:= c;
            fi;
        od;
        unvisited[c]:=1;
        U:= Filtered([1..Y!.nrCells(0)], x -> unvisited[x]=0); 
        Apply(U, x -> distance[x]);
        U:=Minimum(U);
        U:=Positions(distance, U);
        U:=Filtered(U, x -> unvisited[x]=0);
        c:= U[1];
    od;

#Converts the above list into an edge path together with a list of integers (1 or 2) indicating the starting and ending vertex in each edge 
    path:=[vT];
    c:=vT;
    while not c = vS do
        Add(path, previous[c]);
        c:=previous[c];
    od;
    path:=Reversed(path);
    properpath:=[];
    for l in [1..Length(path)-1] do
        Add(properpath, [path[l], path[l+1]]);
    od;
    for m in [1..Length(properpath)] do
        if properpath[m] in strictskeleton then
            properpath[m]:=[Position(skeleton, Set(properpath[m])), 1];
        else
            properpath[m]:=[Position(skeleton, Set(properpath[m])), 2];
        fi;
    od;
    return [distance[vT], properpath];
               
end;
##########################################################
##########################################################
