function [CNR] = CalculateCNR(Signal, TissueModel, Tissue1, Tissue2)
    Rows = size(TissueModel,1);
    Cols = size(TissueModel,2);
    
    s1 = [];
    s2 = [];
    n = [];
    for r = 1:Rows
        for c = 1:Cols
            if strcmp(Tissue1.Name, TissueModel(r,c).Name)
                s1 = [s1 abs(Signal(r,c))];
            elseif strcmp(Tissue2.Name, TissueModel(r,c).Name)
                s2 = [s2 abs(Signal(r,c))];
            elseif strcmp('Empty', TissueModel(r,c).Name)
                n = [n abs(Signal(r,c))];
            end
        end
    end
    
    mean1 = mean(s1);
    mean2 = mean(s2);
    sigmaN = std(n);
    
    CNR = abs(mean1-mean2)/sigmaN;
end
