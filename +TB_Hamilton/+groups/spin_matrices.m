function varargout  = spin_matrices(spin,direction)
% SPIN_MATRICES is used to construct the spin matrix
% Sx = spin_matrices(spin,1); Sy = spin_matrices(spin,2);
% Sz = spin_matrices(spin,3); S0 = spin_matrices(spin,0);
% Ref: https://en.wikipedia.org/wiki/Spin_(physics)
arguments
    spin            (1,1)    double;
    direction     (1,:)    {mustBeInRange(direction,0,3)};
end
d = 2 * spin + 1;
assert(d-fix(d)==0,'The spin must be half-integer');
nargoutchk(length(direction),length(direction));

e = sqrt( 2*(spin+1)*(1:d-1)-(1:d-1).*(2:d) )/2;
for j1 = 1:length(direction)
    varargout{j1} = spin_component(e,d,direction(j1));
end

end


function S = spin_component(e,d,num)
if num == 0
    S = eye(d);
elseif num == 3
    S = diag(d-1:-2:-d+1)/2;
elseif num == 2
    S = full(1i*spdiags([-[0,e]' [e,0]'],[1,-1],d,d));
else
    S =  full(spdiags([[0,e]' [e,0]'],[1,-1],d,d));
end

end