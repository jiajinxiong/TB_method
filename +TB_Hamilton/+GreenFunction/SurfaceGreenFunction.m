function gr = SurfaceGreenFunction(H0,H1,E)
    % H1 : lead --->  scatter region
    E = E+1i*1e-7;
    H1 = H1';
    num = length(H0);
    A = zeros(2*num);   B = eye(2*num);
    A(1:num,num+1:end) = speye(num);
    A(num+1:end,1:num) = -H1';
    A(num+1:end,num+1:end) = E*speye(num)-H0;
    B(num+1:end,num+1:end) = H1;
    [T,S,Q,Z] = qz(A,B);
    lambda = ordeig(T,S);
    ind = zeros(1,2*num);
    ind(abs(lambda)<1) = 1;
    [~,~,~,ZS] = ordqz(T,S,Q,Z,ind);
    ZS = sparse(ZS);    I = speye(num);
    gr = (E * I-H0-H1*ZS(num+1:end,1:num)/ZS(1:num,1:num))\I;
    % gr = (E * speye(num)-H0-ZS(1:num,1:num)\(H1*ZS(num+1:end,1:num)))\speye(num);
end