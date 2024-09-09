%% following the paper PRB 82, 045122 (2010)
clc;clear


% Construct the Hamiltonian

S = TB_Hamilton.groups.spin_matrices(1/2) * 2;
C3z = TB_Hamilton.groups.rotation(2*pi/3,[0,0,1],U=expm(-1i*kron(S(:,:,3),eye(2))*pi/3) );
C2x = TB_Hamilton.groups.rotation(pi,[1,0,0],U=1i*kron(S(:,:,1),S(:,:,3)));
I = TB_Hamilton.groups.rotation(0,[0,0,1],"inversion",true,U = kron(eye(2),S(:,:,3)));
T = TB_Hamilton.groups.time_reversal(3,U=1i*kron(S(:,:,2),eye(2)));
G = TB_Hamilton.groups.generate_group([C3z,C2x,I]);

EC = TB_Hamilton.groups.equivalence_class(G);
EC = EC.keys();
GCT = TB_Hamilton.groups.Character_Table(G);


syms kx ky kz a4 a1 a2 a3 real;

% Gamma Matrix
Gamma1 = kron(S(:,:,1),S(:,:,1)); Gamma2 = kron(S(:,:,2),S(:,:,1)); 
Gamma3 = kron(S(:,:,3),S(:,:,1)); Gamma4 = kron(eye(2),S(:,:,2)); 
Gamma5 = kron(eye(2),S(:,:,3));
Gamma12 = (Gamma1*Gamma2-Gamma2*Gamma1)/(2i);
Gamma13 = (Gamma1*Gamma3-Gamma3*Gamma1)/(2i);
Gamma14 = (Gamma1*Gamma4-Gamma4*Gamma1)/(2i);
Gamma15 = (Gamma1*Gamma5-Gamma5*Gamma1)/(2i);
Gamma23 = (Gamma2*Gamma3-Gamma3*Gamma2)/(2i);
Gamma24 = (Gamma2*Gamma4-Gamma4*Gamma2)/(2i);
Gamma25 = (Gamma2*Gamma5-Gamma5*Gamma2)/(2i);
Gamma34 = (Gamma3*Gamma4-Gamma4*Gamma3)/(2i);
Gamma35 = (Gamma3*Gamma5-Gamma5*Gamma3)/(2i);
Gamma45 = (Gamma4*Gamma5-Gamma5*Gamma4)/(2i);

Tb = GCT.Equiv_Class;
% i = 9; j=10;
% chi = Tb(i+1,:).*Tb(j+1,:);
% TB_Hamilton.groups.Irr_decompose(GCT,chi)

% fun = (randn(1) * kx^3 + randn(1)*ky^3+randn(1)*kz^3+3*randn(1)*kx*ky*kz+...
%     randn(1) * kx^2 * ky + randn(1)*kx^2*kz+randn(1)*ky^2*kz+randn(1)*ky^2*kx+randn(1)*kz^2*kx+randn(1)*kz^2*ky)*eye(4);  %10;
% fun =kz*(kx^2-ky^2)*eye(4);
fun = a1*kx*Gamma1+a2*ky*Gamma1+a3*kx*Gamma2+a4*ky*Gamma2;



simplify(TB_Hamilton.groups.Project_Irr(G,1,fun,[kx,ky,kz],"round_off",5)-fun,'Step',100) % 3kxkykz kz(kx^2-ky^2)
% simplify(TB_Hamilton.groups.Project_Irr(G,12,fun,[kx,ky,kz],"round_off",5)-fun)
% for j1 = 1:length(G)
%     simplify(G(j1).apply(fun,sym_var=[kx,ky,kz])-fun)
% end

% %% all the combinations of the first order terms
% H1 = A1 * (Gamma1 * ky -  Gamma2 * kx) + A0 * 1i*(Gamma15 * ky - Gamma25 * kx) + B1 * Gamma4 * kz + 1i * B0 * Gamma45 * kz;
% H1_ = TB_Hamilton.groups.Project_Irr(G,1,H1,[kx,ky,kz],"round_off",5);
% simplify(H1-H1_)



% %% all the combinations of the second order terms
% H2 = round(randn(1) * Gamma5 + randn(1)*eye(4),2)*kz^2 +round(randn(1) * Gamma5 + randn(1)*eye(4),2)*(kx^2+ky^2) +...
%     A0*(1i  * (2 * kx * ky) * Gamma24 - 1i * (kx^2-ky^2) * Gamma14) + A1*(1i * (2 * kx * ky) * Gamma13+ 1i * (kx^2-ky^2) * Gamma23);
% H2_ = TB_Hamilton.groups.Project_Irr(G,1,H2,[kx,ky,kz],"round_off",5);
% simplify(H2-H2_)


% C3z = TB_Hamilton.groups.rotation(2*pi/3,[0,0,1]);