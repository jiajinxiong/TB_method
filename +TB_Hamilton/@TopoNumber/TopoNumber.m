classdef TopoNumber
    methods(Static)
        function chern = ChernByWilsonLoop(mesh,syst,bandOrEnergy,method)
            arguments
                mesh                pde.FEMesh;
                syst                  TB_Hamilton.Builder;
                bandOrEnergy                 (1,1)   double;
                method   (1,:) char {mustBeMember(method,{'SigleBand','MultiBands','Energy'})} = 'MultiBands';
            end
            points = mesh.Nodes';       ind = mesh.Elements';
            ind_len = length(ind);      chern = 0;
            Hamilton = syst.Hamilton;
            
            for j1 = 1:ind_len
                kxs = points(ind(j1,:),1);    kys = points(ind(j1,:),2);
                Hn1 = full(Hamilton(kxs(1),kys(1)));
                Hn2 = full(Hamilton(kxs(2),kys(2)));
                Hn3 = full(Hamilton(kxs(3),kys(3)));
                [un1,E1] = eig(Hn1); [un2,E2] = eig(Hn2); [un3,E3] = eig(Hn3);
                if method == "SigleBand"
                    Psi_n = det(un1(:,bandOrEnergy)'*un2(:,bandOrEnergy))*det(un2(:,bandOrEnergy)'*un3(:,bandOrEnergy))*...
                        det(un3(:,bandOrEnergy)'*un1(:,bandOrEnergy));
                elseif method == "MultiBands"
                    Psi_n = det(un1(:,1:bandOrEnergy)'*un2(:,1:bandOrEnergy))*det(un2(:,1:bandOrEnergy)'*un3(:,1:bandOrEnergy))*...
                        det(un3(:,1:bandOrEnergy)'*un1(:,1:bandOrEnergy));
                else
                    energy = bandOrEnergy;
                    num1 = length(find(diag(E1) < energy));
                    num2 = length(find(diag(E2) < energy));
                    num3 = length(find(diag(E3) < energy));
                    % if max([num1,num2,num3]) ~= min([num1,num2,num3])
                    %     warning('The number of energy bands below the Fermi surface varies for the three points.');
                    % end
                    num = max([num1,num2,num3]);
                    Psi_n = det(un1(:,1:num)'*un2(:,1:num))*det(un2(:,1:num)'*un3(:,1:num))*...
                        det(un3(:,1:num)'*un1(:,1:num));
                end
                chern = chern -1i*log(Psi_n/abs(Psi_n));
            end
            chern = chern/(2*pi);
        end

        function C = Curvature(syst,kk,num,options)
            arguments
                syst    TB_Hamilton.Builder;
                kk      (2,:)   double;
                num     double;
                options.num_loop        double = 10;
                options.ismulband       double = 0;
            end
            Ham = syst.Hamilton;    nk = size(kk,2);
            Delta = 1e-2;       num_loop = options.num_loop;
            C = zeros(1,nk);
            for j1 = 1:nk
                kx = kk(1,j1);  ky = kk(2,j1);
                Wilson_loop = 1;
                k_loop = [kx+linspace(0,Delta,num_loop),(kx+Delta)*ones(1,num_loop),...
                    kx+linspace(Delta,0,num_loop),kx*ones(1,num_loop);...
                    ky*ones(1,num_loop),ky+linspace(0,Delta,num_loop),...
                    (ky+Delta)*ones(1,num_loop),ky+linspace(Delta,0,num_loop)];
                for j2 = 1:length(k_loop)-1
                    H1 = full(Ham(k_loop(1,j2),k_loop(2,j2)));
                    H2 = full(Ham(k_loop(1,j2+1),k_loop(2,j2+1)));
                    [V1,~] = eig(H1);   [V2,~] = eig(H2);
                    if ~options.ismulband
                        Wilson_loop = Wilson_loop * V1(:,num)'*V2(:,num);
                    else
                        Wilson_loop = Wilson_loop * det(V1(:,1:num)'*V2(:,1:num));
                    end
                end
                H1 = full(Ham(k_loop(1,end),k_loop(2,end)));
                H2 = full(Ham(k_loop(1,1),k_loop(2,1)));
                [V1,~] = eig(H1);   [V2,~] = eig(H2);
                if ~options.ismulband
                    Wilson_loop = Wilson_loop * V1(:,num)'*V2(:,num);
                else
                    Wilson_loop = Wilson_loop * det(V1(:,1:num)'*V2(:,1:num));
                end
                C(j1) = 1i*log(Wilson_loop)/Delta^2;
            end
        end
        function Omega_spin = spin_curvature(syst,Spin,kx,ky,FermiEnergyOrBand,methods)
            arguments
                syst    TB_Hamilton.Builder;
                Spin    (:,:) double
                kx      (1,:) double;
                ky      (1,:) double;
                FermiEnergyOrBand (1,1) double;
                methods   (1,:) char {mustBeMember(methods,{'SigleBand','MultiBands','Energy'})} = 'Energy';
            end
            dim = syst.system_graph.numnodes;
            nx = length(kx); ny = length(ky);
            dy = ky(2)-ky(1);  dx = kx(2) - kx(1);

            assert(length(Spin)==dim,'The dimension of Spin must equal Hamiltonian');

            Omega_spin = zeros(nx,ny);
            parfor j1 = 1:nx
                for j2 = 1:ny
                    H = syst.Hamilton(kx(j1),ky(j2));  [V,E] = eig(full(H)); E = diag(E);
                    H_dx = syst.Hamilton(kx(j1)+dx,ky(j2)); H_dy = syst.Hamilton(kx(j1),ky(j2)+dy);
                    vy = (H_dy-H)/dy;  vx = (H_dx-H)/dx;
                    Omega_ = 0;
                    if methods == "Energy"
                        nband = 1:length(find(E<=FermiEnergyOrBand));
                    elseif methods == "MultiBands"
                        nband = 1:FermiEnergyOrBand;
                    else
                        nband = FermiEnergyOrBand;
                    end
                    for i1 = nband
                        for i2 = [1:i1-1,i1+1:dim]
                            if i2~=i1
                                Delta = E(i1)-E(i2);
                                Omega_ = Omega_ + V(:,i1)'*(vy*Spin+Spin*vy)*V(:,i2)*V(:,i2)'*vx*V(:,i1)/Delta^2*dx*dy;
                            end
                        end
                    end
                    Omega_spin(j1,j2) = Omega_;
                end
            end
            Omega_spin = -imag(Omega_spin);
        end




        function V_ = systFBZ(syst)
            % the point contain the (0,0);
            arguments
                syst    TB_Hamilton.Builder;
            end
            a1 = [syst.vec(1,:),0];   a2 = [syst.vec(2,:),0];
            a3 = [0,0,1];
            volume = abs(dot(a1,cross(a2,a3)));
            b1 = 2*pi*cross(a2,a3)/volume;   b2 = 2*pi*cross(a3,a1)/volume;
            b1 = b1(1:2);   b2 = b2(1:2);
            [x,y] = meshgrid(-3:3,-3:3);
            x = x(:);   y = y(:);
            R_Point = double(x * b1 + y * b2);
            DT = delaunayTriangulation(R_Point);
            [V,r] = voronoiDiagram(DT);
            for j1 = 1:length(r)
                V_ = V(r{j1},:);
                if min(V_(:,1))*max(V_(:,1))<0 && min(V_(:,2))*max(V_(:,2))<0 ...
                        &&max(V_(:))<1e3
                    break;
                end
            end
        end
        function mesh = FBZTri(syst,size_grid)
            arguments
                syst        TB_Hamilton.Builder;
                size_grid   (1,1) double=.5;
            end
            FBZ = TB_Hamilton.TopoNumber.systFBZ(syst);
            region=polyshape(FBZ(:,1),FBZ(:,2));
            tr = triangulation(region);
            model = createpde;
            tnodes = tr.Points';
            telements = tr.ConnectivityList';
            geometryFromMesh(model,tnodes,telements);
            mesh=generateMesh(model,"Hmax",size_grid,"GeometricOrder","linear");
        end
    end
end