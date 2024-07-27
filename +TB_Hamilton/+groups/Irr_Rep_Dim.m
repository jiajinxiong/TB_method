function Rep_D = Irr_Rep_Dim(G)
arguments
    G TB_Hamilton.groups.PointGroupElement;
end
nG = length(G); GEC = TB_Hamilton.groups.equivalence_class(G);
% nEC = length(EC);
Rep_D = sort(decomposeToNSquares(nG,length(GEC.values)));
end




function solutions = decomposeToNSquares(n, numSquares)
    solutions = decompose(n, numSquares, []);
end

function [solution, found] = decompose(n, numSquares, currentSolution)
    if numSquares == 0
        if n == 0
            solution = currentSolution;
            found = true;
        else
            solution = [];
            found = false;
        end
        return;
    end
    
    maxSquareRoot = floor(sqrt(n+1-numSquares));
    
    for i = maxSquareRoot:-1:1
        remaining = n - i^2;
        newSolution = [currentSolution, i];
        [solution, found] = decompose(remaining, numSquares - 1, newSolution);
        if found
            return;
        end
    end
    solution = [];
    found = false;
end