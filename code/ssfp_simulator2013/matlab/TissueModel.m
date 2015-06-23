function [Tissue] = TissueModel(Filename, Tissues, FoSD, MaxFreq)
    % TissueModel genearets a Virtual Tissue Model.
    %
    % TissueModel(Filename, Tissue) generates a Virtual Tissue Model. Model
    % is generated based on supplied tissue masked from file name and the
    % given TissueElements in Tissue.
    %
    % TissueModel(Filename, Tissue, FoSD) generates a Virtual Tissue Model. 
    % Model is generated based on supplied tissue masked from file name and 
    % the given TissueElements in Tissue. This model also has TissueElements
    % with a Gaussian Spread Off-Resonance Freuqencies with a Std. Dev. of 
    % FoSD.
    %
    % TissueModel(Filename, Tissue, FoSD, MaxFreq) generates a Virtual 
    % Tissue Model. Model is generated based on supplied tissue masked from
    % file name and the given TissueElements in Tissue. This model also has 
    % TissueElements with a Gaussian Spread of Off-Resonance Freuqencies 
    % with a Std. Dev. of FoSD. This model also has a linear gradient in
    % the y-dir with a maximum frequncy of MaxFreq.
    %
    
    if nargin < 3
       FoSD = 0;
       MaxFreq = 0;
    end
    
    if nargin < 4
        MaxFreq = 0;
    end
    
    % Tissues (T1, T2, ProtonDensity, Off Resonance Frequency (Hz))
    Empty = Tissues.Empty;
    WhiteMatter = Tissues.WhiteMatter;
    Fat = Tissues.Fat;      

    % Tissue Identification Labels
    black = [0 0 0];        % Empty
    white = [255 255 255];  % Tissue 1
    red = [237 28 36];      % Tissue 2
    
    % Read Tissue Image
    t = imread(Filename);
    
    % Create Tissue Models
    Rows = size(t, 1);
    Cols = size(t, 2);
    Tissue = Empty;
    for r = 1:Rows
        for c = 1:Cols
            if sum(squeeze(t(r,c,:)) == red') == 3
                Tissue(r,c) = Fat;
            elseif sum(squeeze(t(r,c,:)) == white') == 3
                Tissue(r,c) = WhiteMatter;
            else
                Tissue(r,c) = Empty;
            end
        end
    end   
    
    % Generate Off-Resonance Frequency Pulse Spread
    for r = 1:Rows
        for c = 1:Cols
            Tissue(r,c).F0 = Tissue(r,c).F0 + FoSD * randn(1, 1);
        end
    end
    
    % Generate Off-Resonace Frequency Gradient in y-dir
%     phi = linspace(0, MaxFreq, Rows);
    phi = linspace(-MaxFreq/2, MaxFreq/2, Rows);
    for r = 1:Rows
        for c = 1:Cols
            Tissue(r,c).F0 = Tissue(r,c).F0 + phi(r);
        end
    end
end
