classdef Polyatomic
    % POLYATOMIC is a Bravais lattice with an arbitrary number of sites in the basis.
    % 
    % Contains `Monatomic` sublattices.  Note that an instance of ``Polyatomic`` is
    % not itself a `~kwant.builder.SiteFamily`, only its sublattices are.
    % 
    % Parameters
    % ----------
    % prim_vecs : 2d array-like of floats
    %     The primitive vectors of the Bravais lattice
    % basis : 2d array-like of floats
    %     The coordinates of the basis sites inside the unit cell.
    % name : string or sequence of strings, optional
    %     The name of the lattice, or a sequence of the names of all the
    %     sublattices.  If the name of the lattice is given, the names of
    %     sublattices are obtained by appending their number to the name of the
    %     lattice.
    % norbs : int or sequence of ints, optional
    %     The number of orbitals per site on the lattice, or a sequence
    %     of the number of orbitals of sites on each of the sublattices.
    % 
    % Raises
    % ------
    % ValueError
    %     If dimensionalities do not match.
    
    
    properties
        prim_vecs;
        sublattice;
        names;
        orbits;
        basises;
        shape;
    end
    
    methods
        function obj = Polyatomic(prim_vecs, basis, name,orbits)
            arguments
                prim_vecs       double
                basis           double
                name            string
                orbits          double = ones(1,length(name));
            end
            obj.prim_vecs = prim_vecs;
            len = length(name);
            obj.names = name;
            obj.sublattice = {};
            obj.orbits = orbits;
            obj.basises = basis;
            for j1 = 1:len
                obj.sublattice{end+1} = ...
                    TB_Hamilton.Monatomic(prim_vecs, obj.names(j1), basis(j1,:),orbits(j1));
            end
        end
        function obj = set.shape(obj,region)
            if ~isa(region,'polyshape')
                error("Only LatInShape can be used when the shape is well defined.");
            end
            obj.shape = region;
            for j1 = 1:length(obj.sublattice)
                obj.sublattice{j1}.shape = obj.shape;
            end
        end
    end
end

