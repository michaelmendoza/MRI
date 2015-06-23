function [Tissue] = TissueModel1D(Tissues, FoSD, MaxFreq)
    % TissueModel genearets a Virtual Tissue Model.
    %
    % TissueModel(Tissue) generates a Virtual Tissue Model. Model
    % is generated based on supplied tissue masked from file name and the
    % given TissueElements in Tissue.
    %
    % TissueModel(Tissue, FoSD) generates a Virtual Tissue Model. 
    % Model is generated based on supplied tissue masked from file name and 
    % the given TissueElements in Tissue. This model also has TissueElements
    % with a Gaussian Spread Off-Resonance Freuqencies with a Std. Dev. of 
    % FoSD.
    %
    % TissueModel(Tissue, FoSD, MaxFreq) generates a Virtual 
    % Tissue Model. Model is generated based on supplied tissue masked from
    % file name and the given TissueElements in Tissue. This model also has 
    % TissueElements with a Gaussian Spread of Off-Resonance Freuqencies 
    % with a Std. Dev. of FoSD. This model also has a linear gradient in
    % the y-dir with a maximum frequncy of MaxFreq.
    %
    
    if nargin < 2
       FoSD = 0;
       MaxFreq = 0;
    end
    
    if nargin < 3
        MaxFreq = 0;
    end
    
    % Tissues (T1, T2, ProtonDensity, Off Resonance Frequency (Hz))
    Empty = Tissues.Empty;
    WhiteMatter = Tissues.WhiteMatter;
    Fat = Tissues.Fat;      
    
    % Create Tissue Mask
    N0 = 20; Nt = 100;
    t = [zeros(1,N0) ones(1,Nt) zeros(1,N0) 2 * ones(1,Nt) zeros(1,N0)];
    N = length(t);
    
    % Create Tissue Models
    Tissue = Empty;
    for n = 1:N
        if t(n) == 2
          	Tissue(n) = Fat;
      	elseif t(n) == 1
        	Tissue(n) = WhiteMatter;
        else
            Tissue(n) = Empty;
        end
    end   
    
    % Generate Off-Resonance Frequency Pulse Spread
    for n = 1:N
    	Tissue(n).F0 = Tissue(n).F0 + FoSD * randn(1, 1);
    end
    
    % Generate Off-Resonace Frequency Gradient
    %phi = linspace(0, MaxFreq, Nt);
    phi = linspace(-MaxFreq/2, MaxFreq/2, Nt);
    for n = 1:Nt
            w = 20+n;  Tissue(w).F0 = Tissue(w).F0 + phi(n);
            f = 140+n; Tissue(f).F0 = Tissue(f).F0 + phi(n);
    end
end