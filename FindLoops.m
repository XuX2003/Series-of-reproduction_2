function [y, z] = FindLoops(x)
global LoopMat Loops_NodeList NodeLoopNameCounter A EdgeNumMat Edges

A = x;
Nodes = length(A);

%sets L equal to the passed in matrix reshaped into a 1 row vector
L = reshape(x, 1, []);
%creates a vector L, by removing all the non-zero elements of L
L (L == 0) = [];
%number of edges is equal to the length of L
Edges = length(L);

%NodeLoopList name counter
NodeLoopNameCounter = 1;

Loops_NodeList = [];

%builds a matrix where the edges are numbered in increasing order
%initialized it to be the same size as the passed in matrix
%really crappy code, needs to be vectorized and re-written
EdgeNumMat = zeros(Nodes, Nodes);
EdgeNum = 1;
for j = 1: Nodes
    for i = 1: Nodes
        if (A(i, j) ~= 0)
            EdgeNumMat(i, j) = EdgeNum;
            EdgeNum = EdgeNum + 1;
        end
    end
end

%initializes LoopMat to empty
LoopMat = [];

%main program call
%for each of the nodes build a tree
for k = 1: Nodes
    BuildTree (k, k:Nodes, k); %#ok<NASGU>
end

%set return value to the LoopMat
y = Loops_NodeList;
z = LoopMat;


%end main program

%recursive function to build the trees and subtrees
%input parameters are
%R is equal to the current node
%DIR is the directory of nodes that we are trying to connect to
%FromVect is a vector of the nodes we passed through to get to the
%current node
function BuildTree(R, DIR, FromVect)
global A Loops_NodeList NodeLoopNameCounter
%for each node in the directory
for j = DIR
    %if the current node in the directory equals the root node and
    %there is an edge between them (A matrix element is non zero)
    if ((j == R) && (A(j, j) ~= 0))
        %then you have found a self-loop.  Call function BuildLoopMat to
        %attach the loop
        BuildLoopMat([R, j])
        NodeName = strcat('Loop', num2str(NodeLoopNameCounter));
        Loops_NodeList.(NodeName) = [R, j];
        NodeLoopNameCounter = NodeLoopNameCounter + 1;
    end

    %if the current node in the directory is not equal to the current
    %node we are building a tree from but there is an edge between them
    %then attach BuildTree(j DIR - {j}) to R
    if ((j ~= R) && (A(j, R) ~= 0))
        % if the current node equals the root node (FromVect(1)) then
        % you have found a loop
        if (j == FromVect(1))
            %call BuildLoopMat to attach the loop
            BuildLoopMat([FromVect, j])
            NodeName = strcat('Loop', num2str(NodeLoopNameCounter));
            Loops_NodeList.(NodeName) = [FromVect,j];
            NodeLoopNameCounter = NodeLoopNameCounter + 1;
        end
        %build DIR_j by removing j from DIR
        DIR_j = DIR;
        DIR_j(DIR_j == j) = [];
        %call BuildTree with j as current node, and DIR_j as the
        %directory and FromVect = [FromVect, j]
        BuildTree(j, DIR_j, [FromVect, j])
    end
end
%end function BuildTree

%funtion that builds the loop matrix
%LoopVect will contain a listing of the nodes that form a specific loop
function BuildLoopMat (LoopVect)
global LoopMat Edges EdgeNumMat

%initialize a column vector equal in length to the number of edges
NewLoop = zeros(Edges, 1);
%starting at the second node in the list through the last node
for i = 2: length(LoopVect)
    %goes to the position in the EdgeNumMat indicated by the current
    %node (LoopVect(i)) and the previous node (LoopVect(i-1)) and then
    %builds a column vector NewLoop by adding a 1 to the row because that
    %edge is being used
    NewLoop(EdgeNumMat(LoopVect(i), LoopVect(i - 1)), 1) = 1;
end
%adds the NewLoop vector onto the growing LoopMat matrix
LoopMat = [LoopMat, NewLoop];

%end function BuildLoopMat
