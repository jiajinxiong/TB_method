function [n,theta] = rotation_to_angle(R)
% ROTATION_TO_ANGLE converts 3D rotation matrix to axis and angle.

arguments
    R   (:,:)  {mustBeReal};
end
% trace( L(i)*L(j) ) ~ delta_{ij}
% trace( L(i)*L(j)*L(k) ) ~ epsilon_{ijk}
% orther = 0
theta = acos( (trace(R)-1)/2);
if theta==0
    n = [];
else
    n = zeros(1,3); L = TB_Hamilton.groups.L_matrices(3);
    for j1 = 1:3
        n(j1) = real(1j * trace(L(:,:,j1) * R));
    end
    absn = norm(n) * sign(sum(n));

    if abs(absn)<1e-7
        % n is zero for 2-fold rotation
        [V,E] = eig(R);
        num = find(abs(diag(E)-1)<1e-6);
        n = (V(:,num)*sign(sum(V(:,num))))';
    else
        n = n/absn;
    end
    theta = atan2(absn, (trace(R)-1));
end
end