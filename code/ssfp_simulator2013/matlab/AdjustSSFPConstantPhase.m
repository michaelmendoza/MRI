function [ img, img2, img3, img4 ] = AdjustSSFPConstantPhase( img1, img2, img3, img4, theta)
    % Init
    s = {img1, img2, img3, img4};
    
    % Adjust Signal Phase
    for n = 1:4
        s{n} = s{n} * exp(-1i*theta(n));
    end
    
    img = s{1};
    img2 = s{2};
    img3 = s{3};
    img4 = s{4};
end
% 
% function [ img, img2, img3, img4 ] = AdjustSSFPConstantPhase( img1, img2, img3, img4, region)
%     % Init
%     rows = length(img1(:,1));
%     cols = length(img1(1,:));
%     s = {img1, img2, img3, img4};
%     
%     % Estimate Mean Phase of Water Region
%     x1 = region(1); y1 = region(2); x2 = region(3); y2 = region(4);
%     N = (y2 - y1 + 1) * (x2 - x1 +1);
%     theta1 = mean(reshape(angle(img1(y1:y2,x1:x2)),1,N));   
%     theta2 = mean(reshape(angle(img2(y1:y2,x1:x2)),1,N));
%     theta3 = mean(reshape(angle(img3(y1:y2,x1:x2)),1,N));
%     theta4 = mean(reshape(angle(img4(y1:y2,x1:x2)),1,N));
%     
%     Mtheta = [theta1, theta2, theta3, theta4];
% %     Mtheta = [1.5708    0.5230    0.8395    2.3021];
%     
%     % Calculate Constant Phase Difference for Water
%     startIndex = 3;
%     for n = 1:4
%         dTheta(n) = Mtheta(n) - Mtheta(startIndex);
%     end
%     
%     % Adjust Signal Phase
%     for n = 1:4
%         s{n} = s{n} * exp(-1i*dTheta(n));
%     end
%     
%     img = s{1};
%     img2 = s{2};
%     img3 = s{3};
%     img4 = s{4};
% end