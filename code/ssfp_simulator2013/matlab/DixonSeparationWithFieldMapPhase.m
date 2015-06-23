function [ water, fat ] = DixonSeparationWithFieldMapPhase( img1, img2, phi )
    % DixonSeparation separates supplied images into a water and fat 
    % component.
    %
    % [water, fat] = DixonSeparation(img1,img2) uses 2-point Dixon to
    % separate water and fat from two images.
    % 
    % [water, fat] = DixonSeparation(img1,img2,phi) uses 3-point Dixon to
    % separate water and fat from two using the phi as the field map field
    % between img1 and img2
    
    if (nargin < 3)
        water = 0.5 * (img1 + img2);
        fat = 0.5 * (img1 - img2);
    else
        water = 0.5*(img1+img2.*exp(-1i*phi));
        fat = 0.5*(img1-img2.*exp(-1i*phi));
    end
end



