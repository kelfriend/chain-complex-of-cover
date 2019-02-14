LoadPackage("HAP");

##########################################################
##########################################################
ShortestPath:=function(Y,w,vS,vT)
    local i, unvisited, distance, previous, Neighbours, c, n, altdist, U, traceback, p, path, pp, ppp;

    #if not IsHapRegularCWComplex(Y) then
    #    Print("The given complex is not a regular CW-complex.\n");
    #    return fail;
    #fi;

    #if not EvaluateProperty(Y, "dimension") < infinity then
    #    Print("The given complex is not finite.\n");
    #    return fail;
    #fi;

    for i in w do
        if not i>=0 then
            Print("The assigned weights must be non-negative.\n");
            break;
            return fail;
        fi;
    od;

    #1skeleton:=Set(Y!.boundaries[2], x->[x[2],x[3]]);
    #Initialisation step (All the 6's below should be Y!.nrCells(0))
    unvisited:=[1..6]; Apply(unvisited, x -> 0);
    distance:=[1..6]; Apply(distance, x -> infinity); distance[vS]:=0;
    previous:=[]; Add(previous,1);

##########################################################
#Outputs the neighbours of a vertex in the 1-skeleton
    Neighbours:=function(a)
        local j, N;
        N:=[];
        for j in Y do
            if a in j then
                Add(N, Filtered(j, x-> not x=a)[1]);
            fi;
        od;
        return N;
    end;
##########################################################
#Obtains shortest distance from vS to vT as well as a list of previously visited vertices in the optimal path  
    c:=vS;
    while not c = vT do
        for n in Neighbours(c) do
            altdist:= distance[c] + w[Position(Y,Set([c,n]))];
            if altdist < distance[n] then
                distance[n]:= altdist;
                previous[n]:=c;
            fi;
        od;
        unvisited[c]:=1;
        U:= Filtered([1..6], x -> unvisited[x]=0); 
        Apply(U, x -> distance[x]);
        U:=Minimum(U);
        U:=Positions(distance, U);
        U:=Filtered(U, x -> unvisited[x]=0);
        c:= U[1];
    od;

#Turns the list obtained above into a path from vT to vS
# !!! Only works if all nodes up to Length(previous) have been visited !!!
    traceback:=[];
    for p in [1..Length(previous)] do
        Add(traceback, [p,previous[p]]);
    od;
    path:=[]; Add(path, traceback[vT]);
    for pp in [1..6] do
        Add(path, traceback[path[pp][2]]);
    od;
    for ppp in [Position(path, [vS,vS])..Length(path)] do
        Unbind(path[ppp]);
    od;
    return [distance[vT], path];

end;
##########################################################
##########################################################
#test data from wikipedia example

Y:=[[1,2], [1,3], [1,6], [2,3], [2,4], [3,4], [3,6], [4,5], [5,6]];
Apply(Y, x -> Set(x));
w:=[7, 9, 14, 10, 15, 11, 2, 6, 9];
vS:=1;
vT:=5;
