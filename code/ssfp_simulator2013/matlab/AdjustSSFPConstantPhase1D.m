function [ img, img2, img3, img4 ] = AdjustSSFPConstantPhase1D( img1, img2, img3, img4, region)

    % Init
    s = {img1, img2, img3, img4};
    
    % Estimate Mean Phase of Water Region
    x1 = region(1); x2 = region(2);
    N = (x2 - x1 +1);
    theta1 = mean(angle(img1(x1:x2))); 
    theta2 = mean(angle(img2(x1:x2)));   
    theta3 = mean(angle(img3(x1:x2)));   
    theta4 = mean(angle(img4(x1:x2))); 
    
    Mtheta = [theta1, theta2, theta3, theta4];
%     Mtheta = [1.5708    0.5230    0.8395    2.3021];
    
    % Calculate Constant Phase Difference for Water
    startIndex = 3;
    for n = 1:4
        dTheta(n) = Mtheta(n) - Mtheta(startIndex);
    end
    
    % Adjust Signal Phase
    for n = 1:4
        s{n} = s{n} * exp(-1i*dTheta(n));
    end
    
    img = s{1};
    img2 = s{2};
    img3 = s{3};
    img4 = s{4};
end