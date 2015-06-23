% Load Real Data

%% ------------------------------------------------------------------------
% View Real Data
if 1
    filenames = {'Data/SSFP_90Tip_11272013/meas_MID325_SSFPdjp_te_5_ps_0_FID1227.dat'
                 'Data/SSFP_90Tip_11272013/meas_MID330_SSFPdjp_te_6_17_ps_90_FID1232.dat'
                 'Data/SSFP_90Tip_11272013/meas_MID326_SSFPdjp_te_5_ps_180_FID1228.dat'
                 'Data/SSFP_90Tip_11272013/meas_MID331_SSFPdjp_te_6_17_ps_270_FID1233.dat'
                 'Data/SSFP_90Tip_11272013/meas_MID332_SSFPdjp_te_3_83_ps_90_FID1234.dat'
                 'Data/SSFP_90Tip_11272013/meas_MID333_SSFPdjp_te_3_83_ps_270_FID1235.dat'
                 };

    data = {};
    for n=1:6
        n
        img = readMeasDataVB15(filenames{n});                % Access datafile
        img = fftshift(fftn(fftshift(img(:,:,:,1))));        % Convert from k-space to image
%         img = img(:,:,32);  % Look at z-slice
%         img = squeeze(img(147,:,:)); % Look at x-slice
        data{n} = img;                              % Save image slice
        size(img)
    end
end

%% ------------------------------------------------------------------------
% View Real Data
if 0
    filenames = {'Data/SSFP_90Tip_11132013/meas_MID181_SSFPdjp_data_set_3_3d_ps0_te3_25_no_grad_FID472.dat'
                 'Data/SSFP_90Tip_11132013/meas_MID184_SSFPdjp_data_set_3_3d_ps90_te4_42_no_grad_FID475.dat'
                 'Data/SSFP_90Tip_11132013/meas_MID185_SSFPdjp_data_set_3_3d_ps180_te5_58_no_grad_FID476.dat'
                 'Data/SSFP_90Tip_11132013/meas_MID186_SSFPdjp_data_set_3_3d_ps270_te6_75_no_grad_FID477.dat'
                 };

    data = {};
    for n=1:4
        n
        img = readMeasDataVB15(filenames{n});                % Access datafile
        img = fftshift(fftn(fftshift(img(:,:,:,1))));        % Convert from k-space to image
        img = img(:,:,32);  % Look at z-slice
%         img = squeeze(img(147,:,:)); % Look at x-slice
        data{n} = img;                              % Save image slice
        size(img)
    end
end

%% ------------------------------------------------------------------------
% View Real Data
if 0
    filenames = {'Data/SSFP_90Tip_11132013/meas_MID173_SSFPdjp_data_set_2_3d_ps0_te3_25_FID464.dat'
                 'Data/SSFP_90Tip_11132013/meas_MID176_SSFPdjp_data_set_2_3d_ps90_te4_42_FID467.dat'
                 'Data/SSFP_90Tip_11132013/meas_MID177_SSFPdjp_data_set_2_3d_ps180_te5_58_FID468.dat'
                 'Data/SSFP_90Tip_11132013/meas_MID178_SSFPdjp_data_set_2_3d_ps270_te6_75_FID469.dat'
                 };

    data = {};
    for n=1:4
        n
        img = readMeasDataVB15(filenames{n});                % Access datafile
        img = fftshift(fftn(fftshift(img(:,:,:,1))));        % Convert from k-space to image
%         img = img(:,:,32);  % Look at z-slice
        img = squeeze(img(144,:,:)); % Look at x-slice -144!!, 146, 148
        data{n} = img;                                       % Save image slice
        size(img)
        
%         figure(); tmp = fftshift(fftn(fftshift(img(:,:,:,1)))); imshow(abs(tmp(:,:,32)),[]);            % Look at z slice
%         figure(); tmp = fftshift(fftn(fftshift(img(:,:,:,1)))); imshow(abs(squeeze(tmp(128,:,:))),[]);  % Look at y slice
    end
end

%% ------------------------------------------------------------------------
% View Real Data
if 0
    filenames = {'Data/SSFP_90Tip_11132013/meas_MID160_SSFPdjp_data_set_1_phase_cycling_0_te_2_66_FID451.dat'
                 'Data/SSFP_90Tip_11132013/meas_MID161_SSFPdjp_data_set_1_phase_cycling_90_te_4_42_FID452.dat'
                 'Data/SSFP_90Tip_11132013/meas_MID162_SSFPdjp_data_set_1_phase_cycling_180_te_5_58_FID453.dat'
                 'Data/SSFP_90Tip_11132013/meas_MID163_SSFPdjp_data_set_1_phase_cycling_270_te_6_75_FID454.dat'
                 };

    data = {};
    for n=1:4
        n
        img = readMeasDataVB15(filenames{n});       % Access datafile
        img = fftshift(fftn(fftshift(img)));        % Convert from k-space to image
        img = img(:,:,1);                           % Take a middle slice
        data{n} = img;                              % Save image slice
        size(img)
    end
end


%% ------------------------------------------------------------------------
% View Real Data
if 0
    filenames = {'Data/SSFP_90Tip_11122013/SSFPdjp_0deg_TE_325_FID430.dat'
                 'Data/SSFP_90Tip_11122013/SSFPdjp_90deg_TE_442_FID431.dat'
                 'Data/SSFP_90Tip_11122013/SSFPdjp_180deg_TE_558_FID432.dat'
                 'Data/SSFP_90Tip_11122013/SSFPdjp_270deg_TE_675_FID433.dat'
        
                 'Data/SSFP_90Tip_11122013/SSFPdjp_0deg_TE_2_7_FID434.dat'
                 'Data/SSFP_90Tip_11122013/SSFPdjp_72deg_TE_387_FID435.dat'
                 'Data/SSFP_90Tip_11122013/SSFPdjp_144deg_TE_504_FID436.dat'
                 'Data/SSFP_90Tip_11122013/SSFPdjp_216deg_TE_621_FID437.dat'
                 'Data/SSFP_90Tip_11122013/SSFPdjp_288deg_TE_738_FID439.dat'
    
                 'Data/SSFP_90Tip_11122013/SSFPdjp_288deg_TE_730_FID438.dat'
                 };

    data = {};
    for n=1:4
        n
        img = readMeasDataVB15(filenames{n});       % Access datafile
        img = fftshift(fftn(fftshift(img)));        % Convert from k-space to image
        img = img(:,:,1);                           % Take a middle slice
        data{n} = img;                              % Save image slice
        size(img)
    end
end

%% ------------------------------------------------------------------------
% View Real Data
if 0
    filenames = {'Data/SSFP/meas_MID228_SSFPKevin_te3_2_pc0_FID11024.dat'
                 'Data/SSFP/meas_MID229_SSFPKevin_te3_2_pc90_FID11025.dat'
                 'Data/SSFP/meas_MID230_SSFPKevin_te3_2_pc180_FID11026.dat'
                 'Data/SSFP/meas_MID231_SSFPKevin_te3_2_pc270_FID11027.dat'
                 'Data/SSFP/meas_MID232_SSFPKevin_te4_4_pc0_FID11028.dat'
                 'Data/SSFP/meas_MID233_SSFPKevin_te4_4_pc90_FID11029.dat'
                 'Data/SSFP/meas_MID234_SSFPKevin_te4_4_pc180_FID11030.dat'
                 'Data/SSFP/meas_MID235_SSFPKevin_te4_4_pc270_FID11031.dat'
                 'Data/SSFP/meas_MID236_SSFPKevin_te5_6_pc0_FID11032.dat'
                 'Data/SSFP/meas_MID237_SSFPKevin_te5_6_pc90_FID11033.dat'
                 'Data/SSFP/meas_MID238_SSFPKevin_te5_6_pc180_FID11034.dat'
                 'Data/SSFP/meas_MID239_SSFPKevin_te5_6_pc270_FID11035.dat'
                 'Data/SSFP/meas_MID240_SSFPKevin_te6_8_pc0_FID11036.dat'
                 'Data/SSFP/meas_MID241_SSFPKevin_te6_8_pc90_FID11037.dat'
                 'Data/SSFP/meas_MID242_SSFPKevin_te6_8_pc180_FID11038.dat'
                 'Data/SSFP/meas_MID243_SSFPKevin_te6_8_pc270_FID11039.dat'

                 'Data/SSFP/meas_MID189_SSFPKevin_te3_8_pc0_FID10985.dat'
                 'Data/SSFP/meas_MID190_SSFPKevin_te3_8_pc90_FID10986.dat'
                 'Data/SSFP/meas_MID191_SSFPKevin_te3_8_pc180_FID10987.dat'
                 'Data/SSFP/meas_MID192_SSFPKevin_te3_8_pc270_FID10988.dat'
                 'Data/SSFP/meas_MID193_SSFPKevin_te5_pc0_FID10989.dat'
                 'Data/SSFP/meas_MID194_SSFPKevin_te5_pc90_FID10990.dat'
                 'Data/SSFP/meas_MID195_SSFPKevin_te5_pc180_FID10991.dat'
                 'Data/SSFP/meas_MID196_SSFPKevin_te5_pc270_FID10992.dat'
                 'Data/SSFP/meas_MID197_SSFPKevin_te6_2_pc0_FID10993.dat'
                 'Data/SSFP/meas_MID198_SSFPKevin_te6_2_pc90_FID10994.dat'
                 'Data/SSFP/meas_MID199_SSFPKevin_te6_2_pc180_FID10995.dat'
                 'Data/SSFP/meas_MID200_SSFPKevin_te6_2_pc270_FID10996.dat'
                 'Data/SSFP/meas_MID201_SSFPKevin_te7_3_pc0_FID10997.dat'
                 'Data/SSFP/meas_MID202_SSFPKevin_te7_3_pc90_FID10998.dat'
                 'Data/SSFP/meas_MID203_SSFPKevin_te7_3_pc180_FID10999.dat'
                 'Data/SSFP/meas_MID204_SSFPKevin_te7_3_pc270_FID11000.dat'

                 'Data/SSFP/meas_MID211_SSFPKevin_te3_8_pc0_FID11007.dat'
                 'Data/SSFP/meas_MID212_SSFPKevin_te3_8_pc90_FID11008.dat'
                 'Data/SSFP/meas_MID213_SSFPKevin_te3_8_pc180_FID11009.dat'
                 'Data/SSFP/meas_MID214_SSFPKevin_te3_8_pc270_FID11010.dat'
                 'Data/SSFP/meas_MID215_SSFPKevin_te5_pc0_FID11011.dat'
                 'Data/SSFP/meas_MID216_SSFPKevin_te5_pc90_FID11012.dat'
                 'Data/SSFP/meas_MID217_SSFPKevin_te5_pc180_FID11013.dat'
                 'Data/SSFP/meas_MID218_SSFPKevin_te5_pc270_FID11014.dat'
                 'Data/SSFP/meas_MID219_SSFPKevin_te6_2_pc0_FID11015.dat'
                 'Data/SSFP/meas_MID220_SSFPKevin_te6_2_pc90_FID11016.dat'
                 'Data/SSFP/meas_MID221_SSFPKevin_te6_2_pc180_FID11017.dat'
                 'Data/SSFP/meas_MID222_SSFPKevin_te6_2_pc270_FID11018.dat'
                 'Data/SSFP/meas_MID223_SSFPKevin_te7_3_pc0_FID11019.dat'
                 'Data/SSFP/meas_MID224_SSFPKevin_te7_3_pc90_FID11020.dat'
                 'Data/SSFP/meas_MID225_SSFPKevin_te7_3_pc180_FID11021.dat'
                 'Data/SSFP/meas_MID226_SSFPKevin_te7_3_pc270_FID11022.dat'
                 };

    index = [5 6 7 8];
    %index = [1 6 11 16];      % different TE and phase-cycle, TE symmetric data
    %index = [1 5 9 13];       % different TE same dPhi, TE symmetric data
    %index = [1 6 11 16] + 16; % different TE and phase-cycle, TE anti-symmetric data
    %index = [1 6 11 16] + 32; % different TE and phase-cycle, TE anti-symmetric data, Linear-Gradient
    %index = [5 6 7 8] + 32;
    
    data = {};  
    for n=1:4
        img = readMeasDataVB15(filenames{index(n)});    % Access datafile
        img = fftshift(fftn(fftshift(img)));            % Convert from k-space to image
        img = img(:,:,7);                              % Take a middle slice
        data{n} = img;                                  % Save image slice
    end
end

%% ------------------------------------------------------------------------
% View Real Data
if 0
    filenames = {'Data/SSFP_TE/SSFPKevin_TE5_PC18.dat'
                 'Data/SSFP_TE/SSFPKevin_TE5_PC36.dat'
                 'Data/SSFP_TE/SSFPKevin_TE5_PC54.dat'
                 'Data/SSFP_TE/SSFPKevin_TE5_PC72.dat'
                 'Data/SSFP_TE/SSFPKevin_TE5_PC90.dat'
                 'Data/SSFP_TE/SSFPKevin_TE5_PC108.dat'
                 'Data/SSFP_TE/SSFPKevin_TE5_PC126.dat'
                 'Data/SSFP_TE/SSFPKevin_TE5_PC144.dat'
                 'Data/SSFP_TE/SSFPKevin_TE5_PC162.dat'
                 'Data/SSFP_TE/SSFPKevin_TE5_PC180.dat'
                 'Data/SSFP_TE/SSFPKevin_TE5_PC198.dat'
                 'Data/SSFP_TE/SSFPKevin_TE5_PC216.dat'
                 'Data/SSFP_TE/SSFPKevin_TE5_PC234.dat'
                 'Data/SSFP_TE/SSFPKevin_TE5_PC252.dat'
                 'Data/SSFP_TE/SSFPKevin_TE5_PC270.dat'
                 'Data/SSFP_TE/SSFPKevin_TE5_PC288.dat'
                 'Data/SSFP_TE/SSFPKevin_TE5_PC306.dat'
                 'Data/SSFP_TE/SSFPKevin_TE5_PC324.dat'
                 'Data/SSFP_TE/SSFPKevin_TE5_PC342.dat'
                 'Data/SSFP_TE/SSFPKevin_TE5_PC360.dat'
                 };

    data = {};
    for n=1:20
        n
        img = readMeasDataVB15(filenames{n});       % Access datafile
        img = fftshift(fftn(fftshift(img)));        % Convert from k-space to image
        img = img(:,:,7);                           % Take a middle slice
        data{n} = img;                              % Save image slice
    end
end
%% ------------------------------------------------------------------------
% View Real Data
if 0
    filenames = {'Data/SSFP_dPhi/SSFPKevin_TE3_8_PC360.dat'
                 'Data/SSFP_dPhi/SSFPKevin_TE4_PC360.dat'
                 'Data/SSFP_dPhi/SSFPKevin_TE4_2_PC360.dat'
                 'Data/SSFP_dPhi/SSFPKevin_TE4_4_PC360.dat'
                 'Data/SSFP_dPhi/SSFPKevin_TE4_6_PC360.dat'
                 'Data/SSFP_dPhi/SSFPKevin_TE4_8_PC360.dat'
                 'Data/SSFP_dPhi/SSFPKevin_TE5_PC360.dat'
                 'Data/SSFP_dPhi/SSFPKevin_TE5_2_PC360.dat'
                 'Data/SSFP_dPhi/SSFPKevin_TE5_4_PC360.dat'
                 'Data/SSFP_dPhi/SSFPKevin_TE5_6_PC360.dat'
                 'Data/SSFP_dPhi/SSFPKevin_TE5_8_PC360.dat'
                 'Data/SSFP_dPhi/SSFPKevin_TE6_PC360.dat'
                 'Data/SSFP_dPhi/SSFPKevin_TE6_2_PC360.dat'
                 'Data/SSFP_dPhi/SSFPKevin_TE6_4_PC360.dat'
                 'Data/SSFP_dPhi/SSFPKevin_TE6_6_PC360.dat'
                 'Data/SSFP_dPhi/SSFPKevin_TE6_8_PC360.dat'
                 'Data/SSFP_dPhi/SSFPKevin_TE7_PC360.dat'
                 'Data/SSFP_dPhi/SSFPKevin_TE7_2_PC360.dat'
                 'Data/SSFP_dPhi/SSFPKevin_TE7_4_PC360.dat'
                 };

    data = {};
    for n=1:19
        n
        img = readMeasDataVB15(filenames{n});       % Access datafile
        img = fftshift(fftn(fftshift(img)));        % Convert from k-space to image
        img = img(:,:,7);                           % Take a middle slice
        data{n} = img;                              % Save image slice
    end
end
%% ------------------------------------------------------------------------
% View Real Data
if 0
    filenames = {'Data/GRE_Phantom/GRE_TE_5.dat'
                 'Data/GRE_Phantom/GRE_TE_6_17.dat'
                 'Data/GRE_Phantom/GRE_TE_7_34.dat'
                 };

    data = {};
    for n=1:3
        n
        img = readMeasDataVB15(filenames{n});       % Access datafile
        img = fftshift(fftn(fftshift(img)));        % Convert from k-space to image
        img = img(:,:,2);                           % Take a middle slice
        data{n} = img;                              % Save image slice
    end
end
