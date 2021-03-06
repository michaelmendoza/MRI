function [Tissues] = TissueElements()
    % TissueElements() generates a struct of tissue elements. 
    % 
    % Tissue elements are stuct which contain the following properties:
    %   Name - name of type of TissueElement
    %   ProtonDensity - proton density 
    %   T1 - longitudual relaxation time constant
    %   T2 - transverse relaxation time constant
    %   f0 - off-resonance precession frequency

    % Tissues (T1, T2, ProtonDensity, Off Resonance Frequency (Hz))
    Empty       = struct('Name', 'Empty', ...
                         'T1', 0, ...            
                         'T2', 0, ...            
                         'ProtonDensity', 0, ... 
                         'F0', 0);             
    WhiteMatter = struct('Name', 'WhiteMatter', ...
                         'T1', 790/1000, ...         
                         'T2', 92/1000, ...          
                         'ProtonDensity', 2, ... 
                         'F0', 0);             
    Fat         = struct('Name', 'Fat', ...
                         'T1', 270/1000, ...          
                         'T2', 85/1000, ...            
                         'ProtonDensity', 1, ... 
                         'F0', -428);               % -428Hz @ 3 T    
                     
    Tissues = struct('Empty', Empty, 'WhiteMatter', WhiteMatter, 'Fat', Fat');                 
end