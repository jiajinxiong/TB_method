function result = ToSC(V)
    result = zeros(2*size(V));
    result(1:size(V,1),1:size(V,2)) = V;
    result(size(V,1)+1:end,size(V,2)+1:end) = -conj(V);
end

