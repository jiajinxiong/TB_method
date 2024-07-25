function result = mul_(g1,g2)
persistent cache;
if isempty(cache)
    cache = containers.Map('KeyType', 'char', 'ValueType', 'any');
end

% 将输入转换为字符串键
key = [g1,g2];
key = mat2str(key.keyhash);
% 检查缓存中是否已有结果
if isKey(cache, key)
    result = cache(key);
else
    % 计算结果并存储到缓存中
    result = mul(g1,g2);
    cache(key) = result;

    % 控制缓存大小
    maxsize = 1000;
    if cache.Count > maxsize
        % 移除最早添加的键值对
        keys = cache.keys;
        remove(cache, keys{1});
    end
end
end

function g = mul(g1,g2)
[R1,c1,a1,U1] = deal(g1.R,g1.conjugate,g1.antisymmetry,g1.U);
[R2,c2,a2,U2] = deal(g2.R,g2.conjugate,g2.antisymmetry,g2.U);
if isempty(U1) || isempty(U2)
    g_U = [];
elseif c1
    g_U = U1 * conj(U2);
else
    g_U = U1 * U2;
end
g_R = R1 * R2;
g = TB_Hamilton.groups.PointGroupElement(g_R,conjugate = bitxor(c1,c2),antisymmetry = bitxor(a1,a2),U=g_U);
end