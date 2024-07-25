function R = M_Hash()
persistent A
if isempty(A)
    A = memoize(@keyHash);
    A.CacheSize = 10000;
end
R = A;
end