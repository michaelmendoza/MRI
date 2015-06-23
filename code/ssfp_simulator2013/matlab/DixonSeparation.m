function [ water, fat ] = DixonSeparation( img1, img2, img3 )
    % DixonSeparation separates supplied images into a water and fat 
    % component.
    %
    % [water, fat] = DixonSeparation(img1,img2) uses 2-point Dixon to
    % separate water and fat from two images.
    % 
    % [water, fat] = DixonSeparation(img1,img2,img3) uses 3-point Dixon to
    % separate water and fat from three images.
    
    if (nargin < 3)
        water = 0.5 * (img1 + img2);
        fat = 0.5 * (img1 - img2);
    else
        phi = angle(conj(img1).*img3)/2;
        water = 0.5*(img1+img2.*exp(-1i*phi));
        fat = 0.5*(img1-img2.*exp(-1i*phi));
    end
end



