function H = Ham_Ribbon(obj,Ham_,hopping_cell,ind1,Mk)
    arguments
        obj             TB_Hamilton.Builder;
        Ham_            (:,:)   double;
        hopping_cell            double;
        ind1            (:,:)   double;
        Mk.kx           (1,1)   ;
        Mk.ky                    = [];
    end
    ind0 = obj.cell_ind;
    H = Ham_(ind0,ind0);
    for j1 = 1 : size(ind1,2)
        H1 = Ham_(ind0,ind1(:,j1));
        H = H + H1 .* exp(1i * hopping_cell(j1,:) * [Mk.kx;Mk.ky]) ...
            + H1' * exp(-1i * hopping_cell(j1) * [Mk.kx;Mk.ky]);
    end
end

