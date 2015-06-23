% Algorithm Bench with Real Data
% Initialize with Real Data

if 0
    LoadRealData();
end

for n=1:length(data)
%     data{n} = squeeze(data{n}(128,:,:)); % Look at x-slice
    slicedata{n} = squeeze(data{n}(:,:,28));  % Look at z-slice
end

% View Data
%% ------------------------------------------------------------------------
% Plot Magnitude Data
if 0
    figure();
    for n=1:4
        subplot(strcat('14',num2str(n))); imshow(abs(data{n}), []);
    end
end

%% ------------------------------------------------------------------------
% Plot Phase Data
if 0
    figure();
    for n=1:4
        subplot(strcat('14',num2str(n))); imshow(angle(data{n}), []);
    end   
end

% Process Data

%% ------------------------------------------------------------------------
% Take Complex Sum of Four Images and Plot
if 0
    rows = length(data{1}(:,1));
    cols = length(data{1}(1,:));
    avgImg = zeros(rows,cols);
    for r = 1:rows
        for c = 1:cols
            avgImg(r,c) = data{1}(r,c) + data{2}(r,c) + data{3}(r,c) + data{4}(r,c);
        end
    end
    figure();
    subplot('121'); imshow(abs(avgImg),[]);
    subplot('122'); imshow(angle(avgImg),[]);
    title('Complex Sum');  
end

%% ------------------------------------------------------------------------
% Use Dixon Separation Algorithm and Plot
if 0
    [ water, fat ] = DixonSeparation(data{1}, data{2});
    figure();
    subplot('141'); imshow(abs(water),[]);
    subplot('142'); imshow(angle(water),[]);
    subplot('143'); imshow(abs(fat),[]);
    subplot('144'); imshow(angle(fat),[]);
    title('Dixon Separation');   
end

%% ------------------------------------------------------------------------
% Use SSFP Separation Algorithm and Plot
if 1
     theta = -[0 1/4 2/4 3/4] * 2 * pi;
    [ img, img2, img3, img4 ] = AdjustSSFPConstantPhase(slicedata{1}, slicedata{2}, slicedata{3}, slicedata{4}, theta);
    [ water, fat ] = SSFPSeparation(img, img2, img3, img4, 1 );
    
    figure();
    subplot('121'); imshow(abs(water),[]);
    subplot('122'); imshow(abs(fat),[]);
    
%     figure();
%     subplot('221'); imshow(abs(water),[]);
%     subplot('222'); imshow(angle(water),[]);
%     subplot('223'); imshow(abs(fat),[]);
%     subplot('224'); imshow(angle(fat),[]);
%     
%     figure();
%     subplot(1,6,1); imshow(abs(img),[]);
%     subplot(1,6,2); imshow(abs(img2),[]);
%     subplot(1,6,3); imshow(abs(img3),[]);
%     subplot(1,6,4); imshow(abs(img4),[]);
%     subplot(1,6,5); imshow(abs(water),[]);
%     subplot(1,6,6); imshow(abs(fat),[]);
end

%% ------------------------------------------------------------------------
% Use SSFP Separation with Field Map Algorithm and Plot
if 1
     theta = -[0 1/4 2/4 3/4 1/4 3/4] * 2 * pi;
    
    % Adjust Signal Phase
    s = slicedata;
    for n = 1:6
        s{n} = s{n} * exp(-1i*theta(n));
    end
    img = s{1}; img2 = s{2}; img3 = s{3}; img4 = s{4}; img6 = s{5}; img8 = s{5};
    
    [water, fat ] = SSFPSeparationWithFieldMap(img, img2, img3, img4, [], img6, [], img8, 1);
    
    figure();
    subplot('121'); imshow(abs(water),[]);
    subplot('122'); imshow(abs(fat),[]);
    
%     figure();
%     subplot('221'); imshow(abs(water),[]);
%     subplot('222'); imshow(angle(water),[]);
%     subplot('223'); imshow(abs(fat),[]);
%     subplot('224'); imshow(angle(fat),[]);
%     
%     figure();
%     subplot(1,6,1); imshow(abs(img),[]);
%     subplot(1,6,2); imshow(abs(img2),[]);
%     subplot(1,6,3); imshow(abs(img3),[]);
%     subplot(1,6,4); imshow(abs(img4),[]);
%     subplot(1,6,5); imshow(abs(water),[]);
%     subplot(1,6,6); imshow(abs(fat),[]);
end

