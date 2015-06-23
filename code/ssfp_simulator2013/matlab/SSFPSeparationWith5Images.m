function [ water, fat ] = SSFPSeparationWith5Images( img1, img2, img3, img4, img5, dixonmode )
    % SSFPSeparation separates supplied images into water and fat
    % components.
    %
    % [water,fat] = SSFPSeparation(img1,img2,img3,img4,region) uses a 
    % SSFPDixon separation method to separate four images into water and 
    % fat components. 
    %
    % img1   - First image
    % img2   - Second image
    % img3   - Third image
    % img4   - Fourth image
    % region - a boundary box which contains water in all images

    % Init
    rows = length(img1(:,1));
    cols = length(img1(1,:));
    s = {img1, img2, img3, img4, img5};
        
    % Select Top Signal Magnitude
    [top_s, indices] = SelectTop3Pixels(s);
    
    % Reorder images for Dixon Separation (In, Out, In or Out, In, Out)
    %ps = CalculatePhaseState(indices);      % Calculate phase state i.e. In or Out
    %s = OrderImagesForSeparation(s, ps);    % Reorder according to phase state
    
    % Average Signal Magnitude
    for n = 1:3
        mag{n} = abs(top_s{n});
        phase{n} = angle(top_s{n});
    end
    for r = 1:rows
        for c = 1:cols
            norm = mag{1}(r,c) * mag{2}(r,c) * mag{3}(r,c) / 3;
            top_s{1}(r,c) = norm * exp(1i*phase{1}(r,c));
            top_s{2}(r,c) = norm * exp(1i*phase{2}(r,c));
            top_s{3}(r,c) = norm * exp(1i*phase{3}(r,c));
        end
    end
    
%     if dixonmode == 0
%     [ water, fat ] = DixonSeparation( s{1}, s{2} );
%     else
%     [ water, fat ] = DixonSeparation( s{1}, s{2}, s{3} );
%     end
    [water, fat] = DixonSeparationInOrder( s{1}, s{2}, s{3}, s{4}, s{5}, indices{1}, indices{2}, indices{3} );
    
end

function [SelectPixels, SelectIndices] = SelectTop3Pixels(Images)

    % Initialize
    N = length(Images);
    [Ny, Nx] = size(Images{1});
    SelectPixels  = {zeros(Ny,Nx), zeros(Ny,Nx), zeros(Ny,Nx)};
    SelectIndices = {zeros(Ny,Nx), zeros(Ny,Nx), zeros(Ny,Nx)};
    
    % Select top 3 pixels from Phase Cycled Images
    for r = 1:Ny
        for c = 1:Nx
            
            % Select top 3 pixels
            signal = zeros(1, N);
            for n = 1:length(Images)
                signal(n) = Images{n}(r,c);
            end
            
            [mag, indices] = InsertionSort(abs(signal));

            % Save results
            for n = 1:3
                SelectPixels{n}(r,c) = signal(indices(N+1-n));
                SelectIndices{n}(r,c) = indices(N+1-n);
            end
            
        end
    end
end

function [PhaseState] = CalculatePhaseState(Indices)
    N = length(Indices);
    [Ny, Nx] = size(Indices{1});
    for n = 1:N
        for r = 1:Ny
            for c =1:Nx
                Indices{n}(r,c) = mod(Indices{n}(r,c), 2);
            end
        end
    end
    PhaseState = Indices;
end

function [Array, Indices] = InsertionSort(Array)
 
    Indices = 1:numel(Array);
    for n = (2:numel(Array))
 
        value = Array(n);
        index = Indices(n);
        m = n - 1;
 
        while (m >= 1) && (Array(m) > value)
            Array(m+1) = Array(m);
            Indices(m+1) = Indices(m);
            m = m-1;
        end
 
        Array(m+1) = value;
        Indices(m+1) = index;
    end 
end 

function Images = OrderImagesForSeparation(Images, PhaseStates)

    [Ny,Nx] = size(Images{1});
    for r = 1:Ny
        for c = 1:Nx
            
            PhaseCount = PhaseStates{1}(r,c) + PhaseStates{2}(r,c) + PhaseStates{3}(r,c);
            if PhaseCount == 2  % There are 2 In-Phase Images
                
                if PhaseStates{1}(r,c)
                    Pixel1 = Images{1}(r,c);
                    if PhaseStates{2}(r,c)
                       Pixel2 = Images{3}(r,c);
                       Pixel3 = Images{2}(r,c);
                    else
                       Pixel2 = Images{2}(r,c);
                       Pixel3 = Images{3}(r,c);
                    end
                    
                else
                    Pixel1 = Images{2}(r,c);
                    Pixel2 = Images{1}(r,c);
                    Pixel3 = Images{3}(r,c);
                end
                
            else              % There are 2 Out-Phase Images
                
                if PhaseStates{1}(r,c)
                    Pixel1 = Images{2}(r,c);
                    Pixel2 = Images{1}(r,c);
                    Pixel3 = Images{3}(r,c);
                else
                    
                    Pixel1 = Images{1}(r,c);
                    if PhaseStates{2}(r,c)
                       Pixel2 = Images{2}(r,c);
                       Pixel3 = Images{3}(r,c);
                    else
                       Pixel2 = Images{3}(r,c);
                       Pixel3 = Images{2}(r,c);
                    end
                    
                end

            end
            
            % Set newly reorderd pixels
            Images{1}(r,c) = Pixel1;
            Images{2}(r,c) = Pixel2;
            Images{3}(r,c) = Pixel3;
        end
    end
    
end

function [Water, Fat] = ThreePointDixonWithSort(Image1, Image2, Image3, PhaseState1, PhaseState2, PhaseState3)

    [Ny,Nx] = size(Image1);
    Water = zeros(Ny, Nx);
    Fat = zeros(Ny, Nx);
    for r = 1:Ny
        for c = 1:Nx
            
            PhaseCount = PhaseState1(r,c) + PhaseState2(r,c) + PhaseState3(r,c);
            if PhaseCount == 2  % There are 2 In-Phase Images
                
                if PhaseState1(r,c)
                    Pixel1 = Image1(r,c);
                    if PhaseState2(r,c)
                       Pixel2 = Image3(r,c);
                       Pixel3 = Image2(r,c);
                    else
                       Pixel2 = Image2(r,c);
                       Pixel3 = Image3(r,c);
                    end
                    
                else
                    Pixel1 = Image2(r,c);
                    Pixel2 = Image1(r,c);
                    Pixel3 = Image3(r,c);
                end
                
            else              % There are 2 Out-Phase Images
                
                if PhaseState1(r,c)
                    Pixel1 = Image2(r,c);
                    Pixel2 = Image1(r,c);
                    Pixel3 = Image3(r,c);
                else
                    
                    Pixel1 = Image1(r,c);
                    if PhaseState2(r,c)
                       Pixel2 = Image2(r,c);
                       Pixel3 = Image3(r,c);
                    else
                       Pixel2 = Image3(r,c);
                       Pixel3 = Image2(r,c);
                    end
                    
                end

            end
            
            phi = angle(conj(Pixel1)*Pixel3)/2;
            Water(r,c) = 0.5*(Pixel1+Pixel2*exp(-1i*phi));
            Fat(r,c) = 0.5*(Pixel1-Pixel2*exp(-1i*phi));
        end
    end
    

end

function [water, fat] = DixonSeparationInOrder( Image1, Image2, Image3, Image4, Image5, Index1, Index2, Index3 )
    N = 3;
    [Ny, Nx] = size(Image1);
    water = zeros(Ny, Nx);
    fat = zeros(Ny, Nx);
    for n = 1:N
        for r = 1:Ny
            for c =1:Nx
                % Get Index Positions into Index Array
                Index = [0 0 0 0 0];
                Index = Index + (Index1(r,c) == 1) * [1 0 0 0 0] + (Index1(r,c) == 2) * [0 1 0 0 0] + (Index1(r,c) == 3) * [0 0 1 0 0] + (Index1(r,c) == 4) * [0 0 0 1 0] + (Index1(r,c) == 5) * [0 0 0 0 1];
                Index = Index + (Index2(r,c) == 1) * [1 0 0 0 0] + (Index2(r,c) == 2) * [0 1 0 0 0] + (Index2(r,c) == 3) * [0 0 1 0 0] + (Index2(r,c) == 4) * [0 0 0 1 0] + (Index2(r,c) == 5) * [0 0 0 0 1];
                Index = Index + (Index3(r,c) == 1) * [1 0 0 0 0] + (Index3(r,c) == 2) * [0 1 0 0 0] + (Index3(r,c) == 3) * [0 0 1 0 0] + (Index3(r,c) == 4) * [0 0 0 1 0] + (Index3(r,c) == 5) * [0 0 0 0 1];
                    
                if( Index(1) && Index(2) && Index(3) )
                    if( (abs(Image1(r,c)) + abs(Image2(r,c))) > (abs(Image2(r,c)) + abs(Image3(r,c))) )
                        [water(r,c), fat(r,c)] = DixonSeparation(Image1(r,c),Image2(r,c));
                    else
                        [water(r,c), fat(r,c)] = DixonSeparation(Image2(r,c),Image3(r,c));    
                    end
                elseif( Index(1) && Index(2) && Index(4) )
                    [water(r,c), fat(r,c)] = DixonSeparation(Image1(r,c),Image2(r,c));
                elseif( Index(1) && Index(2) && Index(5) )
                    [water(r,c), fat(r,c)] = DixonSeparation(Image1(r,c),Image2(r,c));
                elseif( Index(1) && Index(3) && Index(4) )
                    [water(r,c), fat(r,c)] = DixonSeparation(Image3(r,c),Image4(r,c));
                elseif( Index(1) && Index(3) && Index(5) )
                    water(r,c) = 0;
                    fat(r,c) = 0;
                elseif( Index(3) && Index(4) && Index(5) )
                    if( (abs(Image3(r,c)) + abs(Image4(r,c))) > (abs(Image4(r,c)) + abs(Image5(r,c))) )
                        [water(r,c), fat(r,c)] = DixonSeparation(Image3(r,c),Image4(r,c));
                    else
                        [water(r,c), fat(r,c)] = DixonSeparation(Image4(r,c),Image5(r,c));    
                    end
                elseif( Index(2) && Index(4) && Index(5) ) 
                    [water(r,c), fat(r,c)] = DixonSeparation(Image4(r,c),Image5(r,c));
                elseif( Index(1) && Index(4) && Index(5) )
                    [water(r,c), fat(r,c)] = DixonSeparation(Image4(r,c),Image5(r,c));
                elseif( Index(2) && Index(3) && Index(4) )
                    if( (abs(Image2(r,c)) + abs(Image3(r,c))) > (abs(Image3(r,c)) + abs(Image4(r,c))) )
                        [water(r,c), fat(r,c)] = DixonSeparation(Image2(r,c),Image3(r,c));
                    else
                        [water(r,c), fat(r,c)] = DixonSeparation(Image3(r,c),Image4(r,c));    
                    end
                elseif( Index(2) && Index(3) && Index(5) )  
                    [water(r,c), fat(r,c)] = DixonSeparation(Image2(r,c),Image3(r,c));
                end
                  
                
            end
        end
    end
end