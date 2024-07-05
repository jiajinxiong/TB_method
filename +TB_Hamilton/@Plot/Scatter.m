function  Scatter(pos,rho,propArgs,options2)
arguments
    pos                     (:,2)   double;
    rho                     (:,1)   double;
    propArgs.?matlab.graphics.chart.primitive.BubbleChart
    options2.limsize            (1,:) double =[1,30];
end
propertyCell = namedargs2cell(propArgs);
pos_ = unique(pos,"rows");
Id = [];
if length(pos_) ~= length(pos)
    for j1 = 1:size(pos_,1)
        id = ismember(pos,pos_(j1,:),"rows");
        Id = [Id,find(id)];
    end
    rho = reshape(sum(rho(Id)),[],1);
else
    pos_ = pos;
end
bubblechart(pos_(:,1),pos_(:,2),rho,propertyCell{:});
bubblesize(options2.limsize);
axis equal;
end

