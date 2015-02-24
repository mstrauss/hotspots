function X = move_fix( obj, X )
% MOVE_FIX updates X (edgelist format) when X is actually outside the edge
% dimension.  It does this by randomly selecting an outgoing edge and
% moving X there.

edge_idx  = X(:,1);
delta = X(:,2);

k = size(X,1);  % number of points to be fixed
M = obj.numedges;      % number of edges

edge = zeros(M,k);
Msize = [M k];
edge( sub2ind(Msize, edge_idx', 1:k) ) = 1;
edge = logical(edge);

% find number of neighboring edges
neigh = zeros(k,1);

[S,E] = obj.graph.findedge( edge_idx );
V1 = obj.graph.Edges.EndNodes(:,1);
V2 = obj.graph.Edges.EndNodes(:,2);
L = obj.L;

% start/start
SS = cell2mat( arrayfun( @(s) V1 == s, S', 'UniformOutput', false ) );
% start/end
SE = cell2mat( arrayfun( @(e) V1 == e, E', 'UniformOutput', false ) );
% end/start
ES = cell2mat( arrayfun( @(s) V2 == s, S', 'UniformOutput', false  ) );
% end/end
EE = cell2mat( arrayfun( @(e) V2 == e, E', 'UniformOutput', false  ) );

% assert pairwise disjunctiveness
assert( all( ~all( SS & SE )));
assert( all( ~all( SS & ES )));
assert( all( ~all( SS & EE )));
assert( all( ~all( SE & ES )));
assert( all( ~all( SE & EE )));
assert( all( ~all( ES & EE )));

% candidates & random neighbors
Scand = SS | ES;  Scand( edge ) = 0;
Ecand = SE | EE;  Ecand( edge ) = 0;
% for each column, select a random entry from the ones
Scand = Scand .* rand(Msize);
Ssel = Scand > 0 & Scand ==repmat(max(Scand),Msize(1),1);
Scrit = delta' < 0;
% border edges
border = sum(Ssel)~=1;
Ssel( sub2ind( Msize, edge_idx(Scrit & border), paren( 1:Msize(2), Scrit & border )' ) ) = 1;
%Ssel( edge(:, Scrit & border), Scrit & border ) = 1
%Ssel( edge(:, Scrit & border) ) = 1;
assert( all( size(Ssel) == Msize ));
assert( all( sum(Ssel(:,Scrit))==1 ) );  % must have a single entry per column
[neigh( Scrit ),~] = find( Ssel(:, Scrit) );

Ecand = Ecand .* rand(size(Ecand));
Esel = Ecand > 0 & Ecand==repmat(max(Ecand),size(Ecand,1),1);
Ecrit = (delta > L(edge_idx))';
% border edges
border = sum(Esel)~=1;
Esel( sub2ind( Msize, edge_idx(Ecrit & border), paren( 1:Msize(2), Ecrit & border )' ) ) = 1;
%Esel( edge(:, Ecrit & border), Ecrit & border ) = 1;
assert( all( sum(Esel(:,Ecrit))==1 ) );  % must have a single entry per column
[neigh( Ecrit ),~] = find( Esel(:, Ecrit) );

% displacements
D = zeros(k, 1);
idx = sub2ind( Msize, neigh',1:k);
cond = Ssel(idx) & SS( idx );
D(cond) = -delta(cond);  % delta < 0
cond = Ssel(idx) & ES( idx );
D(cond) = L(neigh(cond)) + delta(cond);  % delta < 0
cond = Esel(idx) & SE( idx );
D(cond) = delta(cond) - L( edge_idx(cond) );
cond = Esel(idx) & EE( idx );
D(cond) = L(neigh(cond)) - ( delta(cond) - L( edge_idx(cond) ) );

% update X
X(:,1) = neigh;
X(:,2) = D;

% fix X recursively if necessary
done = X(:,2) >= 0 & X(:,2) <= L( X(:,1) );
if any(~done)
  X(~done,:) = move_fix( obj, X(~done,:) );
end
end
