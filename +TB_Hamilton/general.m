function lattice = general(prim_vecs, name, basis, orbits)
    % LATTICE create a Bravais lattice of any dimensionality, with any number of sites.
    % Parameters
    % ----------
    % prim_vecs : 2d array-like of floats
    %     The primitive vectors of the Bravais lattice
    % basis : 2d array-like of floats
    %     The coordinates of the basis sites inside the unit cell
    % name : string or sequence of strings
    %     Name of the lattice, or sequence of names of all of the sublattices.
    %     If the name of the lattice is given, the names of sublattices (if any)
    %     are obtained by appending their number to the name of the lattice.
    % norbs : int or sequence of ints, optional
    %     The number of orbitals per site on the lattice, or a sequence
    %     of the number of orbitals of sites on each of the sublattices.
    % 
    % Returns
    % -------
    % lattice : either `Monatomic` or `Polyatomic`
    %     Resulting lattice.
    % 
    % Notes
    % -----
    % This function is largely an alias to the constructors of corresponding
    % lattices.
    arguments
        prim_vecs       double
        name            string
        basis           double = [];
        orbits          double = ones(1,length(name));
    end
    if isempty(basis) & size(prim_vecs,2)==2
        basis = [0,0];
    elseif isempty(basis) & size(prim_vecs,2)==3
        basis = [0,0,0];
    end
    
    if size(basis,1) == 1
        lattice = TB_Hamilton.Monatomic(prim_vecs, name,basis, orbits);
    else
        lattice = TB_Hamilton.Polyatomic(prim_vecs, basis, name, orbits);
    end
end

