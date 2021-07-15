function  [ObjectsThisOrganoid] = FINAL_f_ImageAnalysisPerOperettaOrganoid_cell_count(Label, ch1, ch2, ch3,ch4, ChannelNames, PreviewPath);
%UNTITLED2 Summary of this function goes here
%Script addapted to synculein detection and quantification (excluding TH) 
% different thresholds for synuclein in order to obtain spacificity 

    % vol(ch1, 0, 500) % Alexa 488 >>> filamet synuclein rb 439
    % vol(ch2, 0, 1000) % Alexa 647 >>> TH chk 111
    % vol(ch3, 0, 10000) % HOECHST 33342 >>> Hoechst imtool(max(ch2, [], 3))
    % vol(ch4, 0, 1000) % TRITC >>> total synuclein ms 421
 
    
    %% Initialize variables
    NucleiMask = [];
    TotalSynMask = [];
    FilamntSynMask = [];
    THMask = [];
    
%vol(ch3, 0, 2500)
    ch3BlurSmall = imfilter(double(ch3), fspecial('gaussian', 21, 1), 'symmetric');%vol(ch3BlurSmall)
    ch3BlurBig = imfilter(double(ch3), fspecial('gaussian', 21, 3), 'symmetric');%vol(ch3BlurBig) %%kind of flatfield corraction, to account for different bk in the pic
    ch3DoG = ch3BlurSmall - ch3BlurBig; %vol(ch3DoG, 0, 200, 'hot')
    NucleiMask = ch3DoG > 75; %vol(NucleiMask)%% 50includes all ncleai around the organoids however also includes error that could be taken as nuclei but it is not 
    NucleiMask = bwareaopen(NucleiMask, 20);%vol(NucleiMask)
    ch3LP = imfilter(ch3, fspecial('gaussian', 11, 1), 'symmetric');%vol(ch3LP, 0, 5000, 'hot')
    NucMaskHigh =  (ch3LP > 3500) .* NucleiMask; %vol(NucMaskHigh, 0, 1)
    NucMaskAlive = NucleiMask & ~NucMaskHigh; % vol(NucMaskAlive)
     
   
    %% TH (ch2)
    ch2MedFilt = []; 
    SizeZ = size(ch2, 3);
    parfor p = 1:SizeZ
        ch2MedFilt(:,:,p) = medfilt2(ch2(:,:,p));
    end
  
    THMask = ch2MedFilt > 200; % vol(THMask, 0, 1)
    THDoG = imfilter(ch2, fspecial('gaussian', 11, 1), 'symmetric') - imfilter(ch2, fspecial('gaussian', 31, 10), 'symmetric');
    %vol(THDoG, 0, 300, 'hot')
    THDoGMask = THDoG > 30;
    %vol(THDoGMask)
    THMask = THMask & THDoGMask;
    %vol(THMask, 0, 1)
    %it(max(THMask, [], 3))
    %THMask = THMask & ~NucleiMask; %%%%%% IF THIS IS ON WILL CREATE ERROR IN THE %TH
    THMask = bwareaopen(THMask, 150); %vol(THMask)

    %% TotalSynMask (ch4)  vol(ch4, 0, 1000)
    % different threshold generated to quantify the levels of synuclein
    % based on intensity 
    
    %vol(ch4, 0, 1000)
    ch4MedFilt = [];
    parfor p = 1:size(ch4, 3)
        ch4MedFilt(:,:,p) = medfilt2(ch4(:,:,p));
    end
    %vol(ch4MedFilt, 0, 1500, 'hot')    
    TotalSynMask = ch4MedFilt >  500; % beff 300
    TotalSynMask = bwareaopen(TotalSynMask, 50); %vol(TotalSynMask) 
 
    TotalSynMaskHigh = ch4MedFilt > 1000;
    TotalSynMaskHigh = bwareaopen(TotalSynMaskHigh, 50); %vol(TotalSynMaskHigh) 
 
    TotalSynMaskHigh2 = ch4MedFilt > 1500 ;% 2000;
    TotalSynMaskHigh2 = bwareaopen(TotalSynMaskHigh2, 50); %vol(TotalSynMaskHigh2) 
 
    TotalSynMaskHigh3 = ch4MedFilt > 2000 ;% 2500;
    TotalSynMaskHigh3 = bwareaopen(TotalSynMaskHigh3, 50); %vol(TotalSynMaskHigh3) 
 
    TotalSynMaskHigh4 = ch4MedFilt > 2500 ;% 3000; % too high levels giving 0 as result
    TotalSynMaskHigh4 = bwareaopen(TotalSynMaskHigh4, 50); %vol(TotalSynMaskHigh4) 
 
    %% FilamntSynMask (ch1) vol(ch1, 0, 1000)
    % different thresholds to obtain the specific signal. stronger staining
    % than syncuelin. higher thresholds
    % Ideally synuclein and filament tent to colocalize 
    
    %vol(ch1, 0, 1000)
    ch1MedFilt = [];
    parfor p = 1:size(ch1, 3)
        ch1MedFilt(:,:,p) = medfilt2(ch1(:,:,p));
    end
    
    %vol(ch4MedFilt, 0, 1500, 'hot')    
    FilamntSynMask = ch1MedFilt > 2000; % bef 1000 & 1500
    FilamntSynMask = bwareaopen(FilamntSynMask, 50); %vol(FilamntSynMask) 
    
    FilamntSynMaskHigh = ch1MedFilt > 3000; 
    FilamntSynMaskHigh = bwareaopen(FilamntSynMaskHigh, 50); %vol(FilamntSynMaskHigh) 
    
    FilamntSynMaskHigh2 = ch1MedFilt > 4000; 
    FilamntSynMaskHigh2 = bwareaopen(FilamntSynMaskHigh2, 50); %vol(FilamntSynMaskHigh2) 
 
    FilamntSynMaskHigh3 = ch1MedFilt > 5000; 
    FilamntSynMaskHigh3 = bwareaopen(FilamntSynMaskHigh3, 50); %vol(FilamntSynMaskHigh3) 
 
    FilamntSynMaskHigh4 = ch1MedFilt > 6000; 
    FilamntSynMaskHigh4 = bwareaopen(FilamntSynMaskHigh4, 50); %vol(FilamntSynMaskHigh4) 
    
    FilamntSynMaskHigh5 = ch1MedFilt > 7000; 
    FilamntSynMaskHigh5 = bwareaopen(FilamntSynMaskHigh5, 50); %vol(FilamntSynMaskHigh5) 
    
    %% co-localization and no co-localization syncuelin and TH (not good staining for quantification, but bias in all conditions) 
    
    %co-localization
    DuplexTotalSynTH = TotalSynMask & THMask;
    DuplexTotalSynHighTH = TotalSynMaskHigh & THMask;
    DuplexTotalSynHigh2TH = TotalSynMaskHigh2 & THMask;
    DuplexTotalSynHigh3TH = TotalSynMaskHigh3 & THMask;
    DuplexTotalSynHigh4TH = TotalSynMaskHigh4 & THMask;
   
    % no co-localization 
    DuplexTotalSynNonTH = TotalSynMask & ~THMask;
    DuplexTotalSynHighNonTH = TotalSynMaskHigh & ~THMask;
    DuplexTotalSynHigh2NonTH = TotalSynMaskHigh2 & ~THMask;
    DuplexTotalSynHigh3NonTH = TotalSynMaskHigh3 & ~THMask;
    DuplexTotalSynHigh4NonTH = TotalSynMaskHigh4 & ~THMask;
    
    %% co-localization and no co-localization syncuelin filament and TH (not good staining for quantification, but bias in all conditions) 
    
    %co-localization
    DuplexFilamntSynTH = FilamntSynMask & THMask;
    DuplexFilamntSynHighTH = FilamntSynMaskHigh & THMask;
    DuplexFilamntSynHigh2TH = FilamntSynMaskHigh2 & THMask;
    DuplexFilamntSynHigh3TH = FilamntSynMaskHigh3 & THMask;
    DuplexFilamntSynHigh4TH = FilamntSynMaskHigh4 & THMask;
    DuplexFilamntSynHigh5TH = FilamntSynMaskHigh5 & THMask;
    
    % no co-localization 
    DuplexFilamntSynNonTH = FilamntSynMask & ~THMask;
    DuplexFilamntSynHighNonTH = FilamntSynMaskHigh & ~THMask;
    DuplexFilamntSynHigh2NonTH = FilamntSynMaskHigh2 & ~THMask;
    DuplexFilamntSynHigh3NonTH = FilamntSynMaskHigh3 & ~THMask;
    DuplexFilamntSynHigh4NonTH = FilamntSynMaskHigh4 & ~THMask;
    DuplexFilamntSynHigh5NonTH = FilamntSynMaskHigh5 & ~THMask;
    %% co-localization both synucleins (total and filament) 
    
    DuplexTotalSynFilamntSyn = TotalSynMask & FilamntSynMask;  %vol(DuplexTotalSynFilamntSyn) 
    DuplexTotalSynHighFilamntSynHigh = TotalSynMaskHigh & FilamntSynMaskHigh; %vol(DuplexTotalSynHighFilamntSynHigh) 
  
    %% skipt synuclein notn in TH since Th Ab doesn't properly mark all TH 

    
    %% Previews 
    
    % Scalebar
    imSize = [size(ch1, 1), size(ch1, 2)];
    [BarMask, BarCenter] = f_barMask(200, 0.42, imSize, imSize(1)-200, 200, 25);
    %it(BarMask)
  
    %basic previews 
    PreviewRealSynFilament = cat (3,  imadjust(max(ch1,[],3), [0 ; 0.009],[0 ; 1]));PreviewRealSynFilament = imoverlay2(PreviewRealSynFilament, BarMask, [1 1 1]);
    PreviewRealTH = cat (3,  imadjust(max(ch2,[],3), [0 ; 0.009],[0 ; 1]));    PreviewRealTH = imoverlay2(PreviewRealTH, BarMask, [1 1 1]);
    PreviewRealNuclei = cat (3,  imadjust(max(ch3,[],3), [0 ; 0.009],[0 ; 1]));PreviewRealNuclei = imoverlay2(PreviewRealNuclei, BarMask, [1 1 1]);
    PreviewRealSynucleinTotal = cat (3,  imadjust(max(ch4,[],3), [0 ; 0.009],[0 ; 1]));PreviewRealSynucleinTotal = imoverlay2(PreviewRealSynucleinTotal, BarMask, [1 1 1]);
    % it(PreviewNuclei)

    PreviewTH = imoverlay2(imadjust(max(ch2,[],3),[0 0.05]), bwperim(max(THMask,[],3)), [1 0 0]);
    PreviewTH = imoverlay2(PreviewTH, BarMask, [1 1 1]);
    %imtool(PreviewTH)
    
    PreviewHoechst = imoverlay2(imadjust(max(ch3,[],3),[0 0.075]), bwperim(max(NucleiMask,[],3)), [1 0 0]);
    PreviewHoechst = imoverlay2(PreviewHoechst, BarMask, [1 1 1]);
    % imtool(PreviewHoechst)
    
    PreviewHoechstHigh = imoverlay2(imadjust(max(ch3,[],3),[0 0.075]), bwperim(max(NucMaskHigh,[],3)), [1 0 0]);
    PreviewHoechstHigh = imoverlay2(PreviewHoechstHigh, BarMask, [1 1 1]);
    % imtool(PreviewHoechstHigh)
    
    PreviewHoechstAlive = imoverlay2(imadjust(max(ch3,[],3),[0 0.075]), bwperim(max(NucMaskAlive,[],3)), [1 0 0]);
    PreviewHoechstAlive = imoverlay2(PreviewHoechstAlive, BarMask, [1 1 1]);
    % imtool(PreviewHoechstAlive)
    
    PreviewTotalSynMask = imoverlay2(imadjust(max(ch4, [], 3), [0 0.02]), bwperim(max(TotalSynMask,[],3)), [0 0 1]);
    PreviewTotalSynMask = imoverlay2(PreviewTotalSynMask, BarMask, [1 1 1]);
    %imtool(PreviewTotalSynMask)
    
    PreviewTotalSynMaskHigh = imoverlay2(imadjust(max(ch4, [], 3), [0 0.02]), bwperim(max(TotalSynMaskHigh,[],3)), [0 0 1]);
    PreviewTotalSynMaskHigh = imoverlay2(PreviewTotalSynMaskHigh, BarMask, [1 1 1]);
    %imtool(PreviewTotalSynMask)
    
    PreviewTotalSynMaskHigh2 = imoverlay2(imadjust(max(ch4, [], 3), [0 0.02]), bwperim(max(TotalSynMaskHigh2,[],3)), [0 0 1]);
    PreviewTotalSynMaskHigh2 = imoverlay2(PreviewTotalSynMaskHigh2, BarMask, [1 1 1]);
    %imtool(PreviewTotalSynMask2)
   
    PreviewTotalSynMaskHigh3 = imoverlay2(imadjust(max(ch4, [], 3), [0 0.02]), bwperim(max(TotalSynMaskHigh3,[],3)), [0 0 1]);
    PreviewTotalSynMaskHigh3 = imoverlay2(PreviewTotalSynMaskHigh3, BarMask, [1 1 1]);
    %imtool(PreviewTotalSynMaskHigh3)
   
    PreviewTotalSynMaskHigh4 = imoverlay2(imadjust(max(ch4, [], 3), [0 0.02]), bwperim(max(TotalSynMaskHigh4,[],3)), [0 0 1]);
    PreviewTotalSynMaskHigh4 = imoverlay2(PreviewTotalSynMaskHigh4, BarMask, [1 1 1]);
    %imtool(PreviewTotalSynMaskHigh4)
    
    PreviewFilamntSynMask = imoverlay2(imadjust(max(ch1, [], 3), [0 0.08]), bwperim(max(FilamntSynMask,[],3)), [0 0 1]);
    PreviewFilamntSynMask = imoverlay2(PreviewFilamntSynMask, BarMask, [1 1 1]);
    %imtool(PreviewFilamntSynMask)

    PreviewFilamntSynMaskHigh = imoverlay2(imadjust(max(ch1, [], 3), [0 0.08]), bwperim(max(FilamntSynMaskHigh,[],3)), [0 0 1]);
    PreviewFilamntSynMaskHigh = imoverlay2(PreviewFilamntSynMaskHigh, BarMask, [1 1 1]);
    %imtool(PreviewFilamntSynMask)
    
    PreviewFilamntSynMaskHigh2 = imoverlay2(imadjust(max(ch1, [], 3), [0 0.08]), bwperim(max(FilamntSynMaskHigh2,[],3)), [0 0 1]);
    PreviewFilamntSynMaskHigh2 = imoverlay2(PreviewFilamntSynMaskHigh2, BarMask, [1 1 1]);
    %imtool(PreviewFilamntSynMaskHigh2)
    
    PreviewFilamntSynMaskHigh3 = imoverlay2(imadjust(max(ch1, [], 3), [0 0.08]), bwperim(max(FilamntSynMaskHigh3,[],3)), [0 0 1]);
    PreviewFilamntSynMaskHigh3 = imoverlay2(PreviewFilamntSynMaskHigh3, BarMask, [1 1 1]);
    %imtool(PreviewFilamntSynMaskHigh3)
          
    PreviewFilamntSynMaskHigh4 = imoverlay2(imadjust(max(ch1, [], 3), [0 0.08]), bwperim(max(FilamntSynMaskHigh4,[],3)), [0 0 1]);
    PreviewFilamntSynMaskHigh4 = imoverlay2(PreviewFilamntSynMaskHigh4, BarMask, [1 1 1]);
    %imtool(PreviewFilamntSynMaskHigh4)

       
    PreviewFilamntSynMaskHigh5 = imoverlay2(imadjust(max(ch1, [], 3), [0 0.08]), bwperim(max(FilamntSynMaskHigh5,[],3)), [0 0 1]);
    PreviewFilamntSynMaskHigh5 = imoverlay2(PreviewFilamntSynMaskHigh5, BarMask, [1 1 1]);
    %imtool(PreviewFilamntSynMaskHigh5)

    % previews combined (nuclei + synuclein + Filament)  
    PreviewCombined = cat (3,  imadjust(max(ch4,[],3), [0 ; 0.009],[0 ; 1]) , imadjust(max(ch3,[],3),[0 ; 0.009],[0 ; 1]), imadjust(max(ch1,[],3),[0 ; 0.02],[0 ; 1]));
    % it(PreviewCombined)
  
    %
    IdentityString = [Label.AreaName{:}, '_Idx_', num2str(Label.Idx)];
    imwrite(PreviewTH, [PreviewPath, filesep, IdentityString, '_', 'TH', '.png'])
    
    imwrite(PreviewHoechst, [PreviewPath, filesep, IdentityString, '_', 'Hoechst', '.png'])
    imwrite(PreviewHoechstHigh, [PreviewPath, filesep, IdentityString, '_', 'HoechstHigh', '.png'])
    imwrite(PreviewHoechstAlive, [PreviewPath, filesep, IdentityString, '_', 'HoechstAlive', '.png'])
    
    imwrite(PreviewTotalSynMask, [PreviewPath, filesep, IdentityString, '_', 'TotalSynMask', '.png'])
    imwrite(PreviewTotalSynMaskHigh, [PreviewPath, filesep, IdentityString, '_', 'TotalSynMaskHigh', '.png'])
    imwrite(PreviewTotalSynMaskHigh2, [PreviewPath, filesep, IdentityString, '_', 'TotalSynMaskHigh2', '.png'])
    imwrite(PreviewTotalSynMaskHigh3, [PreviewPath, filesep, IdentityString, '_', 'TotalSynMaskHigh3', '.png'])
    imwrite(PreviewTotalSynMaskHigh4, [PreviewPath, filesep, IdentityString, '_', 'TotalSynMaskHigh4', '.png'])

    imwrite(PreviewFilamntSynMask, [PreviewPath, filesep, IdentityString, '_', 'FilamntSynMask', '.png'])
    imwrite(PreviewFilamntSynMaskHigh, [PreviewPath, filesep, IdentityString, '_', 'FilamntSynMaskHigh', '.png'])
    imwrite(PreviewFilamntSynMaskHigh2, [PreviewPath, filesep, IdentityString, '_', 'FilamntSynMaskHigh2', '.png'])
    imwrite(PreviewFilamntSynMaskHigh3, [PreviewPath, filesep, IdentityString, '_', 'FilamntSynMaskHigh3', '.png'])
    imwrite(PreviewFilamntSynMaskHigh4, [PreviewPath, filesep, IdentityString, '_', 'FilamntSynMaskHigh4', '.png'])
    imwrite(PreviewFilamntSynMaskHigh5, [PreviewPath, filesep, IdentityString, '_', 'FilamntSynMaskHigh5', '.png'])
    
    imwrite(PreviewCombined, [PreviewPath, filesep, IdentityString, '_', 'PreviewCombined', '.png'])
    
    imwrite(PreviewRealSynFilament, [PreviewPath, filesep, IdentityString, '_', 'RealSynFilament', '.png'])
    imwrite(PreviewRealTH, [PreviewPath, filesep, IdentityString, '_', 'RealTH', '.png'])
    imwrite(PreviewRealNuclei, [PreviewPath, filesep, IdentityString, '_', 'RealNuclei', '.png'])
    imwrite(PreviewRealSynucleinTotal, [PreviewPath, filesep, IdentityString, '_', 'RealSynucleinTotal', '.png'])

    %% Feature extraction
    
    % MASK
    % Obtain Masks values and normalised values 
    ObjectsThisOrganoid = table();
    ObjectsThisOrganoid.LabelIdx = {Label.Idx};
    ObjectsThisOrganoid.AreaName = {Label.AreaName};
    %
    ObjectsThisOrganoid.NucMaskSum = sum(NucleiMask(:));
    ObjectsThisOrganoid.NucMaskHighSum = sum(NucMaskHigh(:));
    ObjectsThisOrganoid.NucMaskAliveSum = sum(NucMaskAlive(:));
    % percentage alive/death 
    ObjectsThisOrganoid.PercDeathNuc = sum(NucMaskHigh(:))/sum(NucleiMask(:));
    ObjectsThisOrganoid.PercAliveNuc = sum(NucMaskAlive(:))/sum(NucleiMask(:));
    % control proper quantif nuclei
    ObjectsThisOrganoid.PercNucinNuc = sum(NucleiMask(:))/sum(NucleiMask(:));
    % TH 
    ObjectsThisOrganoid.THMaskSum = sum(THMask(:));
    % Total synuclein in thresholds 
    ObjectsThisOrganoid.TotalSynMaskSum = sum(TotalSynMask(:));
    ObjectsThisOrganoid.TotalSynMaskHighSum = sum(TotalSynMaskHigh(:));
    ObjectsThisOrganoid.TotalSynMaskHigh2Sum = sum(TotalSynMaskHigh2(:));
    ObjectsThisOrganoid.TotalSynMaskHigh3Sum = sum(TotalSynMaskHigh3(:));
    ObjectsThisOrganoid.TotalSynMaskHigh4Sum = sum(TotalSynMaskHigh4(:));
    % Filament synuclein in thresholds 
    ObjectsThisOrganoid.FilamntSynMaskSum = sum(FilamntSynMask(:));
    ObjectsThisOrganoid.FilamntSynMaskHighSum = sum(FilamntSynMaskHigh(:));
    ObjectsThisOrganoid.FilamntSynMaskHigh2Sum = sum(FilamntSynMaskHigh2(:));
    ObjectsThisOrganoid.FilamntSynMaskHigh3Sum = sum(FilamntSynMaskHigh3(:));
    ObjectsThisOrganoid.FilamntSynMaskHigh4Sum = sum(FilamntSynMaskHigh4(:));
    ObjectsThisOrganoid.FilamntSynMaskHigh5Sum = sum(FilamntSynMaskHigh5(:));
    % duplex synuclein & TH 
    ObjectsThisOrganoid.DuplexTotalSynTH = sum(DuplexTotalSynTH(:));
    ObjectsThisOrganoid.DuplexTotalSynHighTH = sum(DuplexTotalSynHighTH(:));
    ObjectsThisOrganoid.DuplexTotalSynHigh2TH = sum(DuplexTotalSynHigh2TH(:));
    ObjectsThisOrganoid.DuplexTotalSynHigh3TH = sum(DuplexTotalSynHigh3TH(:));
    ObjectsThisOrganoid.DuplexTotalSynHigh4TH = sum(DuplexTotalSynHigh4TH(:));
    % duplex synuclein Non TH 
    ObjectsThisOrganoid.DuplexTotalSynNonTH = sum(DuplexTotalSynNonTH(:));
    ObjectsThisOrganoid.DuplexTotalSynHighNonTH = sum(DuplexTotalSynHighNonTH(:));
    ObjectsThisOrganoid.DuplexTotalSynHigh2NonTH = sum(DuplexTotalSynHigh2NonTH(:));
    ObjectsThisOrganoid.DuplexTotalSynHigh3NonTH = sum(DuplexTotalSynHigh3NonTH(:));
    ObjectsThisOrganoid.DuplexTotalSynHigh4NonTH = sum(DuplexTotalSynHigh4NonTH(:));
    % duplex synuclein Filament & TH 
    ObjectsThisOrganoid.DuplexFilamntSynTH = sum(DuplexFilamntSynTH(:));
    ObjectsThisOrganoid.DuplexFilamntSynHighTH = sum(DuplexFilamntSynHighTH(:));
    ObjectsThisOrganoid.DuplexFilamntSynHigh2TH = sum(DuplexFilamntSynHigh2TH(:));
    ObjectsThisOrganoid.DuplexFilamntSynHigh3TH = sum(DuplexFilamntSynHigh3TH(:));
    ObjectsThisOrganoid.DuplexFilamntSynHigh4TH = sum(DuplexFilamntSynHigh4TH(:));
    ObjectsThisOrganoid.DuplexFilamntSynHigh5TH = sum(DuplexFilamntSynHigh5TH(:));
    % duplex synuclein Filament Non TH 
    ObjectsThisOrganoid.DuplexFilamntSynNonTH = sum(DuplexFilamntSynNonTH(:));
    ObjectsThisOrganoid.DuplexFilamntSynHighNonTH = sum(DuplexFilamntSynHighNonTH(:));
    ObjectsThisOrganoid.DuplexFilamntSynHigh2NonTH = sum(DuplexFilamntSynHigh2NonTH(:));
    ObjectsThisOrganoid.DuplexFilamntSynHigh3NonTH = sum(DuplexFilamntSynHigh3NonTH(:));
    ObjectsThisOrganoid.DuplexFilamntSynHigh4NonTH = sum(DuplexFilamntSynHigh4NonTH(:));
    ObjectsThisOrganoid.DuplexFilamntSynHigh5NonTH = sum(DuplexFilamntSynHigh5NonTH(:));
    % duplex both synuclein 
    ObjectsThisOrganoid.DuplexTotalSynFilamntSyn = sum(DuplexTotalSynFilamntSyn(:));
    ObjectsThisOrganoid.DuplexTotalSynHighFilamntSynHigh = sum(DuplexTotalSynHighFilamntSynHigh(:));
    % Averaged to nuclei alive 
    ObjectsThisOrganoid.THMaskSumAlive = sum(THMask(:))/sum(NucMaskAlive(:));
    %
    ObjectsThisOrganoid.TotalSynMaskSumAlive = sum(TotalSynMask(:))/sum(NucMaskAlive(:));
    ObjectsThisOrganoid.TotalSynMaskHighSumAlive = sum(TotalSynMaskHigh(:))/sum(NucMaskAlive(:));
    ObjectsThisOrganoid.TotalSynMaskHigh2SumAlive = sum(TotalSynMaskHigh2(:))/sum(NucMaskAlive(:));
    ObjectsThisOrganoid.TotalSynMaskHigh3SumAlive = sum(TotalSynMaskHigh3(:))/sum(NucMaskAlive(:));
    ObjectsThisOrganoid.TotalSynMaskHigh4SumAlive = sum(TotalSynMaskHigh4(:))/sum(NucMaskAlive(:));
    % 
    ObjectsThisOrganoid.FilamntSynMaskSumAlive = sum(FilamntSynMask(:))/sum(NucMaskAlive(:));
    ObjectsThisOrganoid.FilamntSynMaskHighSumAlive = sum(FilamntSynMaskHigh(:))/sum(NucMaskAlive(:));
    ObjectsThisOrganoid.FilamntSynMaskHigh2SumAlive = sum(FilamntSynMaskHigh2(:))/sum(NucMaskAlive(:));
    ObjectsThisOrganoid.FilamntSynMaskHigh3SumAlive = sum(FilamntSynMaskHigh3(:))/sum(NucMaskAlive(:));
    ObjectsThisOrganoid.FilamntSynMaskHigh4SumAlive = sum(FilamntSynMaskHigh4(:))/sum(NucMaskAlive(:));
    ObjectsThisOrganoid.FilamntSynMaskHigh5SumAlive = sum(FilamntSynMaskHigh5(:))/sum(NucMaskAlive(:));
   
 
    % PIXEL 
            %Get pixel values of the mask  
            NucinNucMask = ch3 .* uint16(NucleiMask);
            NucinNucMaskAlive = ch3 .* uint16(NucMaskAlive);
            THinTHMask = ch2 .* uint16(THMask);
            % synuclein total (pixel mask)
            SyninSynMask = ch4 .* uint16(TotalSynMask);
            SyninSynMaskHigh = ch4 .* uint16(TotalSynMaskHigh);
            SyninSynMaskHigh2 = ch4 .* uint16(TotalSynMaskHigh2);
            SyninSynMaskHigh3 = ch4 .* uint16(TotalSynMaskHigh3);
            SyninSynMaskHigh4 = ch4 .* uint16(TotalSynMaskHigh4);
            % synuclein filament (pixel mask)
            FilamntSyninFilamntSynMask = ch1 .* uint16(FilamntSynMask);
            FilamntSyninFilamntSynMaskHigh = ch1 .* uint16(FilamntSynMaskHigh);
            FilamntSyninFilamntSynMaskHigh2 = ch1 .* uint16(FilamntSynMaskHigh2);
            FilamntSyninFilamntSynMaskHigh3 = ch1 .* uint16(FilamntSynMaskHigh3);
            FilamntSyninFilamntSynMaskHigh4 = ch1 .* uint16(FilamntSynMaskHigh4);
            FilamntSyninFilamntSynMaskHigh5 = ch1 .* uint16(FilamntSynMaskHigh5);
    
    % save sum of pixel value 
    ObjectsThisOrganoid.NucinNucMask = sum(NucinNucMask(:));
    ObjectsThisOrganoid.NucinNucMaskAlive = sum(NucinNucMaskAlive(:));
    ObjectsThisOrganoid.THinTHMask = sum(THinTHMask(:));
    % pixel for synuclein 
    ObjectsThisOrganoid.SyninSynMask = sum(SyninSynMask(:));
    ObjectsThisOrganoid.SyninSynMaskHigh = sum(SyninSynMaskHigh(:));
    ObjectsThisOrganoid.SyninSynMaskHigh2 = sum(SyninSynMaskHigh2(:));
    ObjectsThisOrganoid.SyninSynMaskHigh3 = sum(SyninSynMaskHigh3(:));
    ObjectsThisOrganoid.SyninSynMaskHigh4 = sum(SyninSynMaskHigh4(:));
    % pixel for filament 
    ObjectsThisOrganoid.FilamntSyninFilamntSynMask = sum(FilamntSyninFilamntSynMask(:));
    ObjectsThisOrganoid.FilamntSyninFilamntSynMaskHigh = sum(FilamntSyninFilamntSynMaskHigh(:));
    ObjectsThisOrganoid.FilamntSyninFilamntSynMaskHigh2 = sum(FilamntSyninFilamntSynMaskHigh2(:));
    ObjectsThisOrganoid.FilamntSyninFilamntSynMaskHigh3 = sum(FilamntSyninFilamntSynMaskHigh3(:));
    ObjectsThisOrganoid.FilamntSyninFilamntSynMaskHigh4 = sum(FilamntSyninFilamntSynMaskHigh4(:));
    ObjectsThisOrganoid.FilamntSyninFilamntSynMaskHigh5 = sum(FilamntSyninFilamntSynMaskHigh5(:));
    
    %PIXEL/MASK
    % sum of Pixel value averaged to Mask 
    % AP (averagedPixel)
    ObjectsThisOrganoid.AP_NucinNucMask = sum(NucinNucMask(:)) / sum(NucleiMask(:));
    ObjectsThisOrganoid.AP_NucinNucMaskAlive = sum(NucinNucMaskAlive(:)) / sum(NucMaskAlive(:));
    ObjectsThisOrganoid.AP_THinTHMask = sum(THinTHMask(:))/sum(THMask(:));
    % pixel for synuclein 
    ObjectsThisOrganoid.AP_SyninSynMask = sum(SyninSynMask(:))/sum(TotalSynMask(:));
    ObjectsThisOrganoid.AP_SyninSynMaskHigh = sum(SyninSynMaskHigh(:)) / sum(TotalSynMaskHigh(:));
    ObjectsThisOrganoid.AP_SyninSynMaskHigh2 = sum(SyninSynMaskHigh2(:)) / sum(TotalSynMaskHigh2(:));
    ObjectsThisOrganoid.AP_SyninSynMaskHigh3 = sum(SyninSynMaskHigh3(:)) / sum(TotalSynMaskHigh3(:));
    ObjectsThisOrganoid.AP_SyninSynMaskHigh4 = sum(SyninSynMaskHigh4(:)) / sum(TotalSynMaskHigh4(:));
    % pixel for filament 
    ObjectsThisOrganoid.AP_FilamntSyninFilamntSynMask = sum(FilamntSyninFilamntSynMask(:))/ sum(FilamntSynMask(:));
    ObjectsThisOrganoid.AP_FilamntSyninFilamntSynMaskHigh = sum(FilamntSyninFilamntSynMaskHigh(:))/ sum(FilamntSynMaskHigh(:));
    ObjectsThisOrganoid.AP_FilamntSyninFilamntSynMaskHigh2 = sum(FilamntSyninFilamntSynMaskHigh2(:))/ sum(FilamntSynMaskHigh2(:));
    ObjectsThisOrganoid.AP_FilamntSyninFilamntSynMaskHigh3 = sum(FilamntSyninFilamntSynMaskHigh3(:))/ sum(FilamntSynMaskHigh3(:));
    ObjectsThisOrganoid.AP_FilamntSyninFilamntSynMaskHigh4 = sum(FilamntSyninFilamntSynMaskHigh4(:))/ sum(FilamntSynMaskHigh4(:));
    ObjectsThisOrganoid.AP_FilamntSyninFilamntSynMaskHigh5 = sum(FilamntSyninFilamntSynMaskHigh5(:))/ sum(FilamntSynMaskHigh5(:));
    
    
end

