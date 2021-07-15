function  [ObjectsMAP2isOrganoid] = FINAL_f_ImageAnalysisPerOperettaOrganoid_cell_count(Label, ch1, ch2, ch3,ch4, ChannelNames, PreviewPath);
%UNTITLED2 Summary of MAP2is function goes here
%   Detailed explanation goes here

    % vol(ch1, 0, 500) % Alexa 488 >>> filamet synuclein rb 439
       % vol(ch2, 0, 1000) % Alexa 647 >>> MAP2 chk 111
    % vol(ch3, 0, 10000) % HOECHST 33342 >>> Hoechst imtool(max(ch2, [], 3))
    % vol(ch4, 0, 1000) % TRITC >>> total synuclein ms 421
 
    
    %% Initialize variables
    NucleiMask = [];
    TotalSynMask = [];
    PhospoSynMask = [];
    MAP2Mask = [];
    
    %% Segment nuclei

    ch3BlurSmall = imfilter(double(ch3), fspecial('gaussian', 21, 1), 'symmetric');%vol(ch3BlurSmall)
    ch3BlurBig = imfilter(double(ch3), fspecial('gaussian', 21, 3), 'symmetric');%vol(ch3BlurBig) %%kind of flatfield corraction, to account for different bk in MAP2e pic
    ch3DoG = ch3BlurSmall - ch3BlurBig; %vol(ch3DoG, 0, 200, 'hot')
    NucleiMask = ch3DoG > 75; %vol(NucleiMask)%% 50includes all ncleai around MAP2e organoids however also includes error MAP2at could be taken as nuclei but it is not 
    NucleiMask = bwareaopen(NucleiMask, 20);%vol(NucleiMask)
    ch3LP = imfilter(ch3, fspecial('gaussian', 11, 1), 'symmetric');%vol(ch3LP, 0, 5000, 'hot')
    NucMaskHigh =  (ch3LP > 2500) .* NucleiMask; %vol(NucMaskHigh, 0, 1)
    NucMaskAlive = NucleiMask & ~NucMaskHigh; % vol(NucMaskAlive)
     
   
    %% MAP2 (ch2)
    ch2MedFilt = []; 
    SizeZ = size(ch2, 3);
    parfor p = 1:SizeZ
        ch2MedFilt(:,:,p) = medfilt2(ch2(:,:,p));
    end

    MAP2Mask = ch2MedFilt > 1200; % vol(MAP2Mask, 0, 1)
    MAP2DoG = imfilter(ch2, fspecial('gaussian', 11, 1), 'symmetric') - imfilter(ch2, fspecial('gaussian', 31, 10), 'symmetric');
    %vol(MAP2DoG, 0, 300, 'hot')
    MAP2DoGMask = MAP2DoG > 100;
    %vol(MAP2DoGMask)
    MAP2Mask = MAP2Mask & MAP2DoGMask;
    %vol(MAP2Mask, 0, 1)
    %it(max(MAP2Mask, [], 3))
    MAP2Mask = bwareaopen(MAP2Mask, 200); %vol(MAP2Mask)
    
            
    %% TotalSynMask (ch4)
    
    %vol(ch4, 0, 1000)
    ch4MedFilt = [];
    parfor p = 1:size(ch4, 3)
        ch4MedFilt(:,:,p) = medfilt2(ch4(:,:,p));
    end
    %vol(ch4MedFilt, 0, 1500, 'hot')    
    TotalSynMaskG = ch4MedFilt >  500; % beff 300
    TotalSynMask = bwareaopen(TotalSynMaskG, 10); %vol(TotalSynMask) 
    
    % total synuclein (dots) 
    



    %% PhospoSynMask (ch1)
    
    %vol(ch1, 0, 1000)
    ch1MedFilt = [];
    parfor p = 1:size(ch1, 3)
        ch1MedFilt(:,:,p) = medfilt2(ch1(:,:,p));
    end
    %vol(ch4MedFilt, 0, 1500, 'hot')    
    PhospoSynMaskHigh = ch1MedFilt > 1500; % bef 1000 & 1500
    PhospoSynMaskH = bwareaopen(PhospoSynMaskHigh, 10); %vol(PhospoSynMask) 
%     PhospoSynMaskHextract = bwareaopen(PhospoSynMaskHigh, 500);
    
    PhospoSynMaskLow = ch1MedFilt > 500; % bef 1000 & 1500
    PhospoSynMaskL = bwareaopen(PhospoSynMaskLow, 10); %vol(PhospoSynMask)  
%     PhospoSynMaskLextract = bwareaopen(PhospoSynMaskLow, 500);
%     
%     % extracting dotted filament synuclein 
%     PhospoSynMaskDotsHigh = PhospoSynMaskHigh - PhospoSynMaskHextract;
%     PhospoSynMaskDotsLow = PhospoSynMaskLow - PhospoSynMaskLextract;
%   
    % intracellular dotted synuclein agregatted - co-localization with map2
    %
  
%     i4 = imclose(i4,strel('disk',2));
%     i4a = bwareaopen(i4,(round(pixnum*(2.5/im_size),0))^2); % segment anything larger than 2.5 um square length, then subtract out
%     i4 = i4 - i4a;
   
 %% 
 
     TotalSynMaskPhospoSynMask = PhospoSynMaskHigh & TotalSynMask;
     TotalSynMap2Mask = TotalSynMask & MAP2Mask;
     PhospoSynMap2Mask = PhospoSynMaskHigh & MAP2Mask;
    
%      PhospokDotsLowMap2 = PhospoSynMaskDotsLow & MAP2Mask;
%      PhospokDotsHighMap2 = PhospoSynMaskDotsHigh & MAP2Mask;
%      
     
    %% Previews 
    
    % Scalebar
    imSize = [size(ch1, 1), size(ch1, 2)];
    [BarMask, BarCenter] = f_barMask(200, 0.42, imSize, imSize(1)-200, 200, 25);
    %it(BarMask)
  
    PreviewMAP2 = imoverlay2(imadjust(max(ch2,[],3),[0 0.05]), bwperim(max(MAP2Mask,[],3)), [1 0 0]);
    PreviewMAP2 = imoverlay2(PreviewMAP2, BarMask, [1 1 1]);
    %imtool(PreviewMAP2)
    
    PreviewHoechst = imoverlay2(imadjust(max(ch3,[],3),[0 0.075]), bwperim(max(NucleiMask,[],3)), [1 0 0]);
    PreviewHoechst = imoverlay2(PreviewHoechst, BarMask, [1 1 1]);
    % imtool(PreviewHoechst)
    
    PreviewHoechsHigh = imoverlay2(imadjust(max(ch3,[],3),[0 0.075]), bwperim(max(NucMaskHigh,[],3)), [1 0 0]);
    PreviewHoechsHigh = imoverlay2(PreviewHoechsHigh, BarMask, [1 1 1]);
    % imtool(PreviewHoechsHigh)
    
    PreviewHoechstAlive = imoverlay2(imadjust(max(ch3,[],3),[0 0.075]), bwperim(max(NucMaskAlive,[],3)), [1 0 0]);
    PreviewHoechstAlive = imoverlay2(PreviewHoechstAlive, BarMask, [1 1 1]);
    % imtool(PreviewHoechstAlive)
    
    PreviewTotalSynMask = imoverlay2(imadjust(max(ch4, [], 3), [0 0.02]), bwperim(max(TotalSynMask,[],3)), [0 0 1]);
    PreviewTotalSynMask = imoverlay2(PreviewTotalSynMask, BarMask, [1 1 1]);
    %imtool(PreviewTotalSynMask)
 
    PreviewPhospoSynMaskH = imoverlay2(imadjust(max(ch1, [], 3), [0 0.08]), bwperim(max(PhospoSynMaskH,[],3)), [0 0 1]);
    PreviewPhospoSynMaskH = imoverlay2(PreviewPhospoSynMaskH, BarMask, [1 1 1]);
    %imtool(PreviewPhospoSynMaskH)
    
    PreviewPhospoSynMaskL = imoverlay2(imadjust(max(ch1, [], 3), [0 0.08]), bwperim(max(PhospoSynMaskL,[],3)), [0 0 1]);
    PreviewPhospoSynMaskL = imoverlay2(PreviewPhospoSynMaskL, BarMask, [1 1 1]);
    %imtool(PreviewPhospoSynMaskL)
    
     
  chEmpty = zeros(size(ch1),'uint16');
  
  
%    % only Hoechst
%     PreviewHoechstAlone = cat (3, imadjust(max(chEmpty,[],3), [0 ; 1],[0 ; 1]), imadjust(max(chEmpty,[],3),[0 ; 1],[0 ; 1]), imadjust(max(ch3,[],3),[0 ; 0.2],[0 ; 1]));
%       PreviewHoechstAlone = imoverlay2(PreviewHoechstAlone, BarMask, [1 1 1]);
%       %it(PreviewHoechstAlone)
     
    % only SYN
    PreviewPhospoSynColor = cat (3, imadjust(max(chEmpty,[],3), [0 ; 1],[0 ; 1]), imadjust(max(ch1,[],3),[0 ; 0.02],[0 ; 1]), imadjust(max(chEmpty,[],3),[0 ; 1],[0 ; 1]));
      PreviewPhospoSynColor = imoverlay2(PreviewPhospoSynColor, BarMask, [1 1 1]);
      %it(PreviewPhospoSynColor)
   
       % only SYN
    PreviewTotalSynColor = cat (3, imadjust(max(chEmpty,[],3), [0 ; 1],[0 ; 1]), imadjust(max(ch4,[],3),[0 ; 0.02],[0 ; 1]), imadjust(max(chEmpty,[],3),[0 ; 1],[0 ; 1]));
      PreviewTotalSynColor = imoverlay2(PreviewTotalSynColor, BarMask, [1 1 1]);
      %it(PreviewTotalSynColor)
  
    
    % only MAP2
     PreviewMAP2color = cat (3, imadjust(max(ch2,[],3), [0 ; 0.09],[0 ; 1]), imadjust(max(chEmpty,[],3),[0 ; 1],[0 ; 1]), imadjust(max(chEmpty,[],3),[0 ; 1],[0 ; 1]));
       PreviewMAP2color = imoverlay2(PreviewMAP2color, BarMask, [1 1 1]);
       %it(PreviewMAP2color)
 
 % combined TH MAP2 FOXA2
    PreviewCombinedMAP2SYNSYNP = cat (3, imadjust(max(ch2,[],3), [0 ; 0.1],[0 ; 1]), imadjust(max(ch1,[],3),[0 ; 0.05],[0 ; 1]), imadjust(max(ch4,[],3),[0 ; 0.05],[0 ; 1]));
     PreviewCombinedMAP2SYNSYNP = imoverlay2(PreviewCombinedMAP2SYNSYNP, BarMask, [1 1 1]);
    %it(PreviewCombinedMAP2SYNSYNP)
    

    
    
    IdentityString = [Label.AreaName{:}, '_Idx_', num2str(Label.Idx)];
    imwrite(PreviewMAP2, [PreviewPath, filesep, IdentityString, '_', 'MAP2', '.png'])
    imwrite(PreviewHoechst, [PreviewPath, filesep, IdentityString, '_', 'Hoechst', '.png'])

    imwrite(PreviewHoechsHigh, [PreviewPath, filesep, IdentityString, '_', 'HoechsHigh', '.png'])
    imwrite(PreviewHoechstAlive, [PreviewPath, filesep, IdentityString, '_', 'HoechstAlive', '.png'])
    imwrite(PreviewTotalSynMask, [PreviewPath, filesep, IdentityString, '_', 'TotalSynMask', '.png'])
    
    imwrite(PreviewPhospoSynMaskH, [PreviewPath, filesep, IdentityString, '_', 'PhospoSynMaskH', '.png'])
    imwrite(PreviewPhospoSynMaskL, [PreviewPath, filesep, IdentityString, '_', 'PhospoSynMaskL', '.png'])

     imwrite(PreviewPhospoSynColor, [PreviewPath, filesep, IdentityString, '_', 'PreviewPhospoSynColor', '.png'])
  imwrite(PreviewTotalSynColor, [PreviewPath, filesep, IdentityString, '_', 'PreviewTotalSynColor', '.png'])
  imwrite(PreviewMAP2color, [PreviewPath, filesep, IdentityString, '_', 'PreviewMAP2color', '.png'])

    
    
    
    
    %% Feature extraction
    
    ObjectsMAP2isOrganoid = table();
    ObjectsMAP2isOrganoid.LabelIdx = {Label.Idx};
    ObjectsMAP2isOrganoid.AreaName = {Label.AreaName};
    
    ObjectsMAP2isOrganoid.NucMaskSum = sum(NucleiMask(:));
    ObjectsMAP2isOrganoid.NucMaskHighSum = sum(NucMaskHigh(:));
    ObjectsMAP2isOrganoid.NucMaskAliveSum = sum(NucMaskAlive(:));
    
    ObjectsMAP2isOrganoid.MAP2MaskSum = sum(MAP2Mask(:));
    ObjectsMAP2isOrganoid.TotalSynMaskSum = sum(TotalSynMask(:));
    ObjectsMAP2isOrganoid.PhospoSynMaskSumLow = sum(PhospoSynMaskL(:));
    ObjectsMAP2isOrganoid.PhospoSynMaskSumHigh = sum(PhospoSynMaskH(:));
   
    ObjectsMAP2isOrganoid.MAP2_Phospo = sum(PhospoSynMap2Mask(:));
    ObjectsMAP2isOrganoid.MAP2_TotalSynuclein = sum(TotalSynMap2Mask(:));
    ObjectsMAP2isOrganoid.MAP2_DuplexSynuclein = sum(TotalSynMaskPhospoSynMask(:));
    
       
%     ObjectsMAP2isOrganoid.FilamentDoteHigh = sum(PhospoSynMaskDotsHigh(:));
%     ObjectsMAP2isOrganoid.FilamentDoteLow = sum(PhospoSynMaskDotsLow(:));
%     
%     % intracell 
%     ObjectsMAP2isOrganoid.FilamentDoteHighInMap2 = sum(PhospokDotsHighMap2(:));
%     ObjectsMAP2isOrganoid.FilamentDoteLowInMap2 = sum(PhospokDotsLowMap2(:));
%    

end

