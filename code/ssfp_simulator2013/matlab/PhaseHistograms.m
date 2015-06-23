function PhaseHistograms(img1, img2, img3, img4)
    % PhaseHistograms(img1,img2,img3,img4) plots histogram of the phase of 
    % supplied images

    rows = length(img1(:,1));
    cols = length(img1(1,:));

    img1 = reshape(img1,1,rows*cols);
    img2 = reshape(img2,1,rows*cols);
    img3 = reshape(img3,1,rows*cols);
    img4 = reshape(img4,1,rows*cols);
    
    figure();
    subplot(221);
    hist(angle(img1),1000);
    subplot(222);
    hist(angle(img2),1000);
    subplot(223);
    hist(angle(img3),1000);
    subplot(224);
    hist(angle(img4),1000);
end