function Partition(syst,options)
    arguments
        syst        TB_Hamilton.GreenFunction.Green_Function;
        options.?matlab.graphics.chart.primitive.Scatter;
    end
    D = syst.Distance;
    propertyCell = namedargs2cell(options);
    for j1 = 1:max(D)
        Id = find(D==j1);
        pos = syst.GF_pos(Id,:);
        scatter(pos(:,1),pos(:,2),"filled",propertyCell{:});hold on;
    end
    axis equal;
end

