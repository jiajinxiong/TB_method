function H = Ham_Bulk(obj,Ham_,hopping_cell,ind1,Mk)
    % arguments
    % 
    % end

    ind0 = obj.cell_ind;
    pos0 = obj.ham_pos(ind0,:) * [kx;ky];
    [Pos0,Pos1] = meshgrid(pos0,pos0);
    Pos = Pos0-Pos1;
    H = Ham_(ind0,ind0).* exp(1i * Pos);
    for j1 = 1:size(hopping_cell,1)
        H1 = Ham_(ind0,ind1(:,j1));
        pos1 = obj.ham_pos(ind1(:,j1),:) * [kx;ky];
        [Pos0,Pos1] = meshgrid(pos0,pos1);
        Pos = Pos0-Pos1;
        H = H + H1 .* exp(1i * Pos);
    end
    H = (H+H')/2;
end

