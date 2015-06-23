function PlotMRImage( Image1, Image2, Image3, Image4, Image5, Image6 )
    % PlotMRImage plots the magnitude and phase for an MR image

    if nargin < 2
        PlotMagnitudePhaseImage(Image1);
    elseif nargin < 3
        PlotMagnitiudePhaseFor2Images(Image1, Image2);
    elseif nargin < 4
        PlotMagnitiudePhaseFor3Images(Image1, Image2, Image3);
    elseif nargin < 5
        PlotMagnitiudePhaseFor4Images(Image1, Image2, Image3, Image4);
    elseif nargin < 6
        PlotMagnitiudePhaseFor5Images(Image1, Image2, Image3, Image4, Image5);
    elseif nargin < 7
        PlotMagnitiudePhaseFor6Images(Image1, Image2, Image3, Image4, Image5, Image6);
    end

end

function PlotMagnitudePhaseImage(S)
	figure();
  	subplot(1,2,1);
 	imshow(abs(S),[]);         % Signal Magnitude
    title('Magitude Image');

 	subplot(1,2,2);
 	imshow(angle(S),[]);       % Signal Phase
    title('Phase Image');

end

function PlotMagnitiudePhaseFor2Images(S, S2)
    figure();
    subplot(2,2,1);
 	imshow(abs(S),[]);         % Signal Magnitude
    title('MagitudeImage');
 	subplot(2,2,2);
 	imshow(angle(S),[]);       % Signal Phase
    title('Phase Image');
    subplot(2,2,3);
 	imshow(abs(S2),[]);         % Signal Magnitude
 	subplot(2,2,4);
 	imshow(angle(S2),[]);       % Signal Phase
end

function PlotMagnitiudePhaseFor3Images(S, S2, S3)
    figure();
    subplot(3,2,1);
 	imshow(abs(S),[]);         % Signal Magnitude
    title('Magitude Image');
 	subplot(3,2,2);
 	imshow(angle(S),[]);       % Signal Phase
    title('Phase Image');
    subplot(3,2,3);
 	imshow(abs(S2),[]);         % Signal Magnitude
 	subplot(3,2,4);
 	imshow(angle(S2),[]);       % Signal Phase
    subplot(3,2,5);
 	imshow(abs(S3),[]);         % Signal Magnitude
 	subplot(3,2,6);
 	imshow(angle(S3),[]);       % Signal Phase
end

function PlotMagnitiudePhaseFor4Images(S, S2, S3, S4)
    figure();
    subplot(4,2,1);
 	imshow(abs(S),[]);          % Signal Magnitude
    title('Magitude Image');
 	subplot(4,2,2);
 	imshow(angle(S),[]);        % Signal Phase
    title('Phase Image');
    subplot(4,2,3);
 	imshow(abs(S2),[]);         % Signal Magnitude
 	subplot(4,2,4);
 	imshow(angle(S2),[]);       % Signal Phase
    subplot(4,2,5);
 	imshow(abs(S3),[]);         % Signal Magnitude
 	subplot(4,2,6);
 	imshow(angle(S3),[]);       % Signal Phase
    subplot(4,2,7);
    imshow(abs(S4),[]);         % Signal Magnitude
 	subplot(4,2,8);
 	imshow(angle(S4),[]);       % Signal Phase
end

function PlotMagnitiudePhaseFor5Images(S, S2, S3, S4, S5)
    figure();
    subplot(5,2,1);
 	imshow(abs(S),[]);          % Signal Magnitude
    title('Magitude Image');
 	subplot(5,2,2);
 	imshow(angle(S),[]);        % Signal Phase
    title('Phase Image');
    subplot(5,2,3);
 	imshow(abs(S2),[]);         % Signal Magnitude
 	subplot(5,2,4);
 	imshow(angle(S2),[]);       % Signal Phase
    subplot(5,2,5);
 	imshow(abs(S3),[]);         % Signal Magnitude
 	subplot(5,2,6);
 	imshow(angle(S3),[]);       % Signal Phase
    subplot(5,2,7);
    imshow(abs(S4),[]);         % Signal Magnitude
 	subplot(5,2,8);
 	imshow(angle(S4),[]);       % Signal Phase
    subplot(5,2,9);
 	imshow(abs(S5),[]);         % Signal Magnitude
 	subplot(5,2,10);
 	imshow(angle(S5),[]);       % Signal Phase
end

function PlotMagnitiudePhaseFor6Images(S, S2, S3, S4, S5, S6)
    figure();
    subplot(6,2,1);
 	imshow(abs(S),[]);          % Signal Magnitude
    title('Magitude Image');
 	subplot(6,2,2);
 	imshow(angle(S),[]);        % Signal Phase
    title('Phase Image');
    subplot(6,2,3);
 	imshow(abs(S2),[]);         % Signal Magnitude
 	subplot(6,2,4);
 	imshow(angle(S2),[]);       % Signal Phase
    subplot(6,2,5);
 	imshow(abs(S3),[]);         % Signal Magnitude
 	subplot(6,2,6);
 	imshow(angle(S3),[]);       % Signal Phase
    subplot(6,2,7);
    imshow(abs(S4),[]);         % Signal Magnitude
 	subplot(6,2,8);
 	imshow(angle(S4),[]);       % Signal Phase
    subplot(6,2,9);
 	imshow(abs(S5),[]);         % Signal Magnitude
 	subplot(6,2,10);
 	imshow(angle(S5),[]);       % Signal Phase
    subplot(6,2,11);
 	imshow(abs(S6),[]);         % Signal Magnitude
 	subplot(6,2,12);
 	imshow(angle(S6),[]);       % Signal Phase
end

