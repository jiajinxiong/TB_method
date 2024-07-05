function result = SCPair(hopping)
    arguments   
        hopping     double
    end
    % result = [0 hopping;hopping' 0];
    result = zeros(2*size(hopping));
    result(1:size(hopping,1),size(hopping,2)+1:end) = hopping;
    result(size(hopping,1)+1:end,1:size(hopping,2)) = -hopping';
end

