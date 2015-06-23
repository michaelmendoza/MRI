function [ water, fat ] = SSFPSeparationWithFieldMap( img1, img2, img3, img4, img5, img6, img7, img8, dixonmode )
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
    s = {img1, img2, img3, img4};
        
    % Select Top Signal Magnitude
    [top_s, indices] = SelectTop3Pixels(s);
    
    if dixonmode == 0
        [water, fat] = DixonSeparationInOrder( s{1}, s{2}, s{3}, s{4}, img5, img6, img7, img8, indices{1}, indices{2}, indices{3} );
    elseif dixonmode == 1
        [water, fat] = DixonSeparationInOrderDebug( s{1}, s{2}, s{3}, s{4}, img5, img6, img7, img8, indices{1}, indices{2}, indices{3} );
    else
        [water, fat] = DixonSeparationInOrderDebug2( s{1}, s{2}, s{3}, s{4}, img5, img6, img7, img8, indices{1}, indices{2}, indices{3} );
    end    
    
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

function [water, fat] = DixonSeparationInOrder( Image1, Image2, Image3, Image4, Image5, Image6, Image7, Image8, Index1, Index2, Index3 )
    N = 3;
    [Ny, Nx] = size(Image1);
    water = zeros(Ny, Nx);
    fat = zeros(Ny, Nx);
    for n = 1:N
        for r = 1:Ny
            for c =1:Nx
                % Get Index Positions into Index Array
                Index = [0 0 0 0];
                Index = Index + (Index1(r,c) == 1) * [1 0 0 0] + (Index1(r,c) == 2) * [0 1 0 0] + (Index1(r,c) == 3) * [0 0 1 0] + (Index1(r,c) == 4) * [0 0 0 1];
                Index = Index + (Index2(r,c) == 1) * [1 0 0 0] + (Index2(r,c) == 2) * [0 1 0 0] + (Index2(r,c) == 3) * [0 0 1 0] + (Index2(r,c) == 4) * [0 0 0 1];
                Index = Index + (Index3(r,c) == 1) * [1 0 0 0] + (Index3(r,c) == 2) * [0 1 0 0] + (Index3(r,c) == 3) * [0 0 1 0] + (Index3(r,c) == 4) * [0 0 0 1];
                  
                % Go through combinations    
                if( Index(1) && Index(2) && Index(3) )
                    if( (abs(Image1(r,c)) + abs(Image2(r,c))) > (abs(Image2(r,c)) + abs(Image3(r,c))) )
                        [water(r,c), fat(r,c)] = DixonSeparation(Image1(r,c),Image2(r,c));
                    else
                        [water(r,c), fat(r,c)] = DixonSeparation(Image2(r,c),Image3(r,c));    
                    end
                elseif( Index(1) && Index(2) && Index(4))
                  	[water(r,c), fat(r,c)] = DixonSeparation(Image1(r,c),Image2(r,c));
                elseif( Index(1) && Index(3) && Index(4))
                    [water(r,c), fat(r,c)] = DixonSeparation(Image3(r,c),Image4(r,c));
                elseif( Index(2) && Index(3) && Index(4))
                    if( (abs(Image2(r,c)) + abs(Image3(r,c))) > (abs(Image3(r,c)) + abs(Image4(r,c))) )
                        [water(r,c), fat(r,c)] = DixonSeparation(Image2(r,c),Image3(r,c));
                    else
                        [water(r,c), fat(r,c)] = DixonSeparation(Image3(r,c),Image4(r,c));    
                    end
                else
                   'Error!' 
                end
                
            end
        end
    end
   
    
end

function [water, fat] = DixonSeparationInOrderDebug( Image1, Image2, Image3, Image4, Image5, Image6, Image7, Image8, Index1, Index2, Index3 )
    N = 3;
    [Ny, Nx] = size(Image1);
    water = zeros(Ny, Nx);
    fat = zeros(Ny, Nx);
    top1 = zeros(Ny, Nx);   % First In-Out
    top2 = zeros(Ny, Nx);   % Second Out-In
    top3 = zeros(Ny, Nx);   % left over
    phase = zeros(Ny, Nx);
    for n = 1:N
        for r = 1:Ny
            for c =1:Nx
                % Get Index Positions into Index Array
                Index = [0 0 0 0];
                Index = Index + (Index1(r,c) == 1) * [1 0 0 0] + (Index1(r,c) == 2) * [0 1 0 0] + (Index1(r,c) == 3) * [0 0 1 0] + (Index1(r,c) == 4) * [0 0 0 1];
                Index = Index + (Index2(r,c) == 1) * [1 0 0 0] + (Index2(r,c) == 2) * [0 1 0 0] + (Index2(r,c) == 3) * [0 0 1 0] + (Index2(r,c) == 4) * [0 0 0 1];
                Index = Index + (Index3(r,c) == 1) * [1 0 0 0] + (Index3(r,c) == 2) * [0 1 0 0] + (Index3(r,c) == 3) * [0 0 1 0] + (Index3(r,c) == 4) * [0 0 0 1];
                  
                % Go through combinations    
                if( Index(1) && Index(2) && Index(3) )
                    if( (abs(Image1(r,c)) + abs(Image2(r,c))) > (abs(Image2(r,c)) + abs(Image3(r,c))) )
                        top1(r,c) = Image1(r,c);
                        top2(r,c) = Image2(r,c);
                        top3(r,c) = Image3(r,c);
                        phase(r,c) = angle(conj(Image6(r,c)).*Image2(r,c)) / 2;
                    else
                        top1(r,c) = Image3(r,c);
                        top2(r,c) = Image2(r,c); 
                        top3(r,c) = Image1(r,c);
                        phase(r,c) = angle(conj(Image6(r,c)).*Image2(r,c)) / 2;
                    end
                elseif( Index(1) && Index(2) && Index(4))
                    top1(r,c) = Image1(r,c);
                    top2(r,c) = Image2(r,c);  
                    top3(r,c) = Image4(r,c);
                    phase(r,c) = angle(conj(Image6(r,c)).*Image2(r,c)) / 2;
                elseif( Index(1) && Index(3) && Index(4))
                        top1(r,c) = Image3(r,c);
                        top2(r,c) = Image4(r,c);
                        top3(r,c) = Image1(r,c);
                        phase(r,c) = angle(conj(Image8(r,c)).*Image4(r,c)) / 2;
                elseif( Index(2) && Index(3) && Index(4))
                    if( (abs(Image2(r,c)) + abs(Image3(r,c))) > (abs(Image3(r,c)) + abs(Image4(r,c))) )
                        top1(r,c) = Image3(r,c);
                        top2(r,c) = Image2(r,c); 
                        top3(r,c) = Image4(r,c);
                        phase(r,c) = angle(conj(Image6(r,c)).*Image2(r,c)) / 2;
                    else
                        top1(r,c) = Image3(r,c);
                        top2(r,c) = Image4(r,c); 
                        top3(r,c) = Image2(r,c);
                        phase(r,c) = angle(conj(Image8(r,c)).*Image4(r,c)) / 2;
                    end
                else
                   'Error!' 
                end
                
            end
        end
    end
    
	for r = 1:Ny
        for c = 1:Nx
            norm = abs(top1(r,c)) + abs(top2(r,c)) + abs(top3(r,c)) / 3;
            [water(r,c), fat(r,c)] = DixonSeparationWithFieldMapPhase( norm * exp(1i*angle(top1(r,c))) , norm * exp(1i*angle(top2(r,c))), phase(r,c) );
        end
    end
end

function [water, fat] = DixonSeparationInOrderDebug2( Image1, Image2, Image3, Image4, Image5, Image6, Image7, Image8, Index1, Index2, Index3 )
    N = 3;
    [Ny, Nx] = size(Image1);
    water = zeros(Ny, Nx);
    fat = zeros(Ny, Nx);
    top1 = zeros(Ny, Nx);   % First In-Out
    top2 = zeros(Ny, Nx);   % Second Out-In
    top3 = zeros(Ny, Nx);   % left over
    for n = 1:N
        for r = 1:Ny
            for c =1:Nx
                % Get Index Positions into Index Array
                Index = [0 0 0 0];
                Index = Index + (Index1(r,c) == 1) * [1 0 0 0] + (Index1(r,c) == 2) * [0 1 0 0] + (Index1(r,c) == 3) * [0 0 1 0] + (Index1(r,c) == 4) * [0 0 0 1];
                Index = Index + (Index2(r,c) == 1) * [1 0 0 0] + (Index2(r,c) == 2) * [0 1 0 0] + (Index2(r,c) == 3) * [0 0 1 0] + (Index2(r,c) == 4) * [0 0 0 1];
                Index = Index + (Index3(r,c) == 1) * [1 0 0 0] + (Index3(r,c) == 2) * [0 1 0 0] + (Index3(r,c) == 3) * [0 0 1 0] + (Index3(r,c) == 4) * [0 0 0 1];
                  
                % Go through combinations    
                if( Index(1) && Index(2) && Index(3) )
                    if( (abs(Image1(r,c)) + abs(Image2(r,c))) > (abs(Image2(r,c)) + abs(Image3(r,c))) )
                        top1(r,c) = Image1(r,c);
                        top2(r,c) = Image2(r,c);
                        top3(r,c) = Image3(r,c);
                        phase(r,c) = angle(conj(Image6(r,c)).*Image2(r,c)) / 2;
                    else
                        top1(r,c) = Image3(r,c);
                        top2(r,c) = Image2(r,c); 
                        top3(r,c) = Image1(r,c);
                        phase(r,c) = angle(conj(Image6(r,c)).*Image2(r,c)) / 2;
                    end
                elseif( Index(1) && Index(2) && Index(4))
                    if( (abs(Image1(r,c)) + abs(Image2(r,c))) > (abs(Image1(r,c)) + abs(Image4(r,c))) )
                      	top1(r,c) = Image1(r,c);
                        top2(r,c) = Image2(r,c);  
                        top3(r,c) = Image4(r,c);
                        phase(r,c) = angle(conj(Image6(r,c)).*Image2(r,c)) / 2;
                    else
                      	top1(r,c) = Image1(r,c);
                        top2(r,c) = Image4(r,c);  
                        top3(r,c) = Image2(r,c);
                        phase(r,c) = angle(conj(Image8(r,c)).*Image4(r,c)) / 2;
                    end

                elseif( Index(1) && Index(3) && Index(4))
                    if( (abs(Image1(r,c)) + abs(Image4(r,c))) > (abs(Image3(r,c)) + abs(Image4(r,c))) )
                        top1(r,c) = Image1(r,c);
                        top2(r,c) = Image4(r,c);
                        top3(r,c) = Image3(r,c);  
                        phase(r,c) = angle(conj(Image8(r,c)).*Image4(r,c)) / 2;
                    else
                        top1(r,c) = Image3(r,c);
                        top2(r,c) = Image4(r,c);
                        top3(r,c) = Image1(r,c); 
                        phase(r,c) = angle(conj(Image8(r,c)).*Image4(r,c)) / 2;
                    end

                elseif( Index(2) && Index(3) && Index(4))
                    if( (abs(Image2(r,c)) + abs(Image3(r,c))) > (abs(Image3(r,c)) + abs(Image4(r,c))) )
                        top1(r,c) = Image3(r,c);
                        top2(r,c) = Image2(r,c); 
                        top3(r,c) = Image4(r,c);
                        phase(r,c) = angle(conj(Image6(r,c)).*Image2(r,c)) / 2;
                    else
                        top1(r,c) = Image3(r,c);
                        top2(r,c) = Image4(r,c); 
                        top3(r,c) = Image2(r,c);
                        phase(r,c) = angle(conj(Image8(r,c)).*Image4(r,c)) / 2;
                    end
                else
                   'Error!' 
                end
                
            end
        end
    end
    
	for r = 1:Ny
        for c = 1:Nx
            norm = abs(top1(r,c)) + abs(top2(r,c)) + abs(top3(r,c)) / 3;
            [water(r,c), fat(r,c)] = DixonSeparation( norm * exp(1i*angle(top1(r,c))) , norm * exp(1i*angle(top2(r,c))) );
        end
    end
end





