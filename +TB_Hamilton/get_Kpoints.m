function [K_Lines,K_labels] = get_Kpoints(High_sym_Points,Nk)
    arguments
       High_sym_Points double ;
       Nk (:,1)  double;
    end
    if ~(length(Nk)~=1 || size(High_sym_Points,1)-1-length(Nk)~=0)
       error('The number of k points and internal points are not matched')
    end
    if isscalar(Nk)
       Nk = Nk * ones(size(High_sym_Points,1)-1,1);
    end
    K_labels = zeros(size(High_sym_Points,1),1);
    K_Lines = zeros(sum(Nk),size(High_sym_Points,2));
    for j1 = 1:size(High_sym_Points,1)-1
       for j2 = 1:size(High_sym_Points,2)
          K_Lines(sum(Nk(1:j1-1))+1:sum(Nk(1:j1)),j2) = linspace(High_sym_Points(j1,j2),High_sym_Points(j1+1,j2),Nk(j1));
       end
       K_labels(j1) = sum(Nk(1:j1-1))+1;
    end
    K_labels(end) = sum(Nk);
 end