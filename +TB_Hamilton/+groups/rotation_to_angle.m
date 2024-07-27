function [n,theta] = rotation_to_angle(R)
% ROTATION_TO_ANGLE converts 3D rotation matrix to axis and angle.

arguments
    R   (3,3)  {mustBeReal};
end
    n = zeros(1,3); 
    L = TB_Hamilton.groups.L_matrices(3);
    for j1 = 1:3
        n(j1) = real(1j * trace(L(:,:,j1) * R));
    end
    absn = norm(n) * sign(sum(n));

    if abs(absn)<1e-5
        % n is zero for 2-fold rotation
        [V,E] = eig(R);
        num = find(abs(diag(E)-1)<1e-6,1);
        % if 
        n = V(:,num);
        n(abs(n)<1e-5) = 0;
        n = n * sign(sum(n)+1e-4);
        if abs(sum(n))<1e-4
            n = n * sign(n(find(n,1)));
        end
    else
        n = n/absn;
    end
    theta = atan2(absn, (trace(R)-1));
    if abs(abs(theta)-pi)<1e-4
       theta = abs(theta); 
    end
end
% end