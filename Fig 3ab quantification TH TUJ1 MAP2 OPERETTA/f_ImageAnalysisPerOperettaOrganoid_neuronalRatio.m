function  [ObjectsThisOrganoid] = f_ImageAnalysisPerOperettaOrganoid_cell_count_C(Label, ch1, ch2, ch3,ch4, ChannelNames, PreviewPath);
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    % vol(ch1, 0, 500) % Alexa 488 >>> TH
       % vol(ch2, 0, 1000) % Alexa 647 >>> MAP2
    % vol(ch3, 0, 10000) % HOECHST 33342 >>> Hoechst imtool(max(ch2, [], 3))
    % vol(ch4, 0, 1000) % TRITC >>> TUJ1
 
    
    %% Initialize variables
    NucleiMask = [];
    THMask = [];
    TUJ1Mask = [];
    MAP2Mask = [];
    
    %% Segment nuclei

%vol(ch3, 0, 2500)
    ch3BlurSmall = imfilter(double(ch3), fspecial('gaussian', 21, 1), 'symmetric');%vol(ch3BlurSmall)
    ch3BlurBig = imfilter(double(ch3), fspecial('gaussian', 21, 3), 'symmetric');%vol(ch3BlurBig) %%kind of flatfield corraction, to account for different bk in the pic
    ch3DoG = ch3BlurSmall - ch3BlurBig; %vol(ch3DoG, 0, 200, 'hot')
    NucleiMask = ch3DoG > 80; %vol(NucleiMask)%% 50includes all ncleai around the organoids however also includes error that could be taken as nuclei but it is not 
    NucleiMask = bwareaopen(NucleiMask, 20);%vol(NucleiMask)
    ch3LP = imfilter(ch3, fspecial('gaussian', 11, 1), 'symmetric');%vol(ch3LP, 0, 5000, 'hot')
    NucMaskHigh =  (ch3LP > 2000) .* NucleiMask; %vol(NucMaskHigh, 0, 1)
    NucMaskAlive = NucleiMask & ~NucMaskHigh; % vol(NucMaskAlive)
     
   %% MAP2 (ch2)647
    
    %vol(ch2, 0, 500)
    ch2MedFilt = []; 
    SizeZ = size(ch2, 3);
    parfor p = 1:SizeZ
        ch2MedFilt(:,:,p) = medfilt2(ch2(:,:,p));
    end
    
    %vol(ch3MedFilt, 0, 1000, 'hot')
    MAP2Mask = ch2MedFilt > 2000; % vol(map2Mask, 0, 1) % bef 250
    MAP2DoG = imfilter(ch2, fspecial('gaussian', 11, 1), 'symmetric') - imfilter(ch2, fspecial('gaussian', 31, 10), 'symmetric');
    %vol(map2DoG, 0, 300, 'hot')
    MAP2DoGMask = MAP2DoG > 300; %bef 50
    %vol(map2DoGMask)
    MAP2Mask = MAP2Mask & MAP2DoGMask;
    %vol(map2Mask, 0, 1)
    MAP2Mask = bwareaopen(MAP2Mask, 150); %vol(MAP2Mask)
    
    %% TH (ch1)
    ch1MedFilt = []; 
    SizeZ = size(ch1, 3);
    parfor p = 1:SizeZ
        ch1MedFilt(:,:,p) = medfilt2(ch1(:,:,p));
    end

    THMask = ch1MedFilt > 500; % 300 for 0.5 % vol(THMask, 0, 1)
    THDoG = imfilter(ch1, fspecial('gaussian', 11, 1), 'symmetric') - imfilter(ch1, fspecial('gaussian', 31, 10), 'symmetric');
    %vol(THDoG, 0, 300, 'hot')
    THDoGMask = THDoG > 300;
    %vol(THDoGMask)
    THMask = THMask & THDoGMask;
    %vol(THMask, 0, 1)
    %it(max(THMask, [], 3))
    THMask = bwareaopen(THMask, 150); %vol(THMask)

% %     %% TH skeleton3D EPFL
% %     disp('Start skel')
% %     tic
% %     skelTH = Skeleton3D(THMask);
% %     toc
% %     disp('Skel done')
% %     %vol(skelTH, 0, 1)
% %     [AdjacencyMatrixTH, nodeTH, linkTH] = Skel2Graph3D(skelTH,0);                       
% %     %imtool(AdjacencyMatrixTH, [])
% %     NodeTH = zeros(size(THMask), 'uint8');
% %     NodeIdxs = vertcat(nodeTH(:).idx);
% %     NodeTH(NodeIdxs) = 1;
% %     %vol(NodeTH)    
% %     if size(NodeIdxs, 1) == 0
% %         return
% %     end
% %     NodeTHPreview = uint8(skelTH) + NodeTH + uint8(THMask); 
% %     NodeTHPreview2D = max(NodeTHPreview, [], 3);
% %     %it(NodeTHPreview2D)
% %     %vol(NodeTHPreview, 0, 3, 'jet')    
% %     NodeDegreeVectorTH = sum(AdjacencyMatrixTH, 1);
% % 
% %     ZeroNodeExplanationNeeded = 0;
% %     if ZeroNodeExplanationNeeded
% %         ZeroNodes = find(NodeDegreeVectorTH == 0);
% %         ZeroNodesLinIdx = vertcat(nodeTH(ZeroNodes).idx);
% %         ZeroNodeMask = zeros(size(THMaskClipped), 'uint8');
% %         ZeroNodeMask(ZeroNodesLinIdx) = 1; %vol(ZeroNodeMask)
% %         NodePreviewZeroCase = uint8(skelTH) + NodeMaskTH + 10*uint8(ZeroNodeMask) + uint8(THMask);
% %     end  
% % 
% %     %% TH Fragmentation
% % 
% %     % Define structuring element for surface detection
% %     Conn6 = strel('sphere', 1); % 6 connectivity
% %     % Detect surface
% %     SurfaceTH = THMask & ~(imerode(THMask, Conn6));
% %     %vol(SurfaceTH)
    
    %% TUJ1 (ch4)
    
     ch1MedFilt = []; 
    SizeZ = size(ch4, 3);
    parfor p = 1:SizeZ
        ch1MedFilt(:,:,p) = medfilt2(ch4(:,:,p));
    end
    TUJ1Mask = ch1MedFilt > 1000; % vol(THMask, 0, 1)
    TUJ1TDoG = imfilter(ch4, fspecial('gaussian', 11, 1), 'symmetric') - imfilter(ch4, fspecial('gaussian', 31, 10), 'symmetric');
    %vol(THDoG, 0, 300, 'hot')
    TUJ1TDoGMask = TUJ1TDoG >500;
    %vol(THDoGMask)
    TUJ1Mask = TUJ1Mask & TUJ1TDoGMask;
    %vol(THMask, 0, 1)
    %it(max(THMask, [], 3))
    TUJ1Mask = bwareaopen(TUJ1Mask, 150); %vol(THMask)
    
% %     %% Perinuclear Volume (to detect amount of cells positive for TH)
% %     
% %     %vol(NucleiMask)
% %     NucleiMaskSingleCells = f_RemoveBigObjects (NucleiMask, 10000); 
% %     NucDil = imdilate(imdilate(NucleiMaskSingleCells, strel('disk', 4)), strel('sphere',1));
% %     NucPerim = logical(NucDil) & ~logical(NucleiMaskSingleCells);
% %     %vol(NucPerim)
% %     THMaskinNucPerim = THMask & NucPerim;% vol(THMaskinNucPerim)

    %% TUJ1 Mask in non-TH
    MaskTHTUJ1 =  THMask & TUJ1Mask;
    MaskTUJ1MAP2 =  MAP2Mask  & TUJ1Mask;  
    MaskMAP2TH = MAP2Mask & THMask ;
    
    %vol(TUJ1MaskNonTH,0,1)        
    %vol(TUJ1Perim,0,1)   
    %vol(TUJ1MaskNonTH,0,1)   
% %     %% Percent TH pos
% %     %split perinuc
% %     D = bwdist(NucleiMaskSingleCells);
% %     %vol(D, 0, 20, 'hot')
% %     %it(D(:,:,1))
% %     disp('start watershed')
% %     tic
% %     W = watershed(D);
% %     toc
% %     disp('watershed done')
% %     %vol(W)
% %     NucPerimStencil = uint16(W) .* uint16(imreconstruct(logical(imdilate(NucPerim, strel('disk', 1))), logical(NucleiMaskSingleCells))); % This line was causing the error 20171207 % Function imreconstruct expected MARKER and MASK to have the same class.
% %     %vol(NucPerimStencil)
% %     %vol(NucPerim)
% %     %vol(NucleiMaskSingleCells)
% %     %toto = imreconstruct(logical(NucPerim), logical(NucleiMaskSingleCells));
% %        
% %     PeriNucMask = logical(NucPerimStencil);
% %     PeriNucMask = bwareaopen(PeriNucMask, 500);
% %     %vol(PeriNucMask)
% %     
% %     PerinucLM = bwlabeln(PeriNucMask);%vol(PerinucLM); vol(uint16(PeriNucMask) .* uint16(THMask), 0,1); vol(THMask +2*PeriNucMask)
% %     PeriNucObjects = regionprops('table', PerinucLM, double(THMask), 'PixelValues');
% %     THproportions = rowfun(@(x) sum(x{:})/length(x{:}), PeriNucObjects, 'InputVariables', 'PixelValues');
% %     THPos = array2table(table2array(THproportions) > 0.01);
% %     THPos.Properties.VariableNames(end) = {'THpos'};
% %     PeriNucObjects = [PeriNucObjects, THproportions, THPos];
% %     PeriNucObjects.Properties.VariableNames(end-1) = {'THproportion'};%{'PixelValues', 'THproportion', 'THpos'};
% %     PeriNucObjectsCompact = PeriNucObjects(:, {'THproportion','THpos'});
% %     THPercent = (sum(PeriNucObjectsCompact.THpos)/height(PeriNucObjectsCompact))*100;
% %     
     
    %% Previews 
    
    % Scalebar
    imSize = [size(ch1, 1), size(ch1, 2)];
    [BarMask, BarCenter] = f_barMask(200, 0.42, imSize, imSize(1)-200, 200, 25);
    %it(BarMask)
    
    PreviewHoechst = imoverlay2(imadjust(max(ch3,[],3),[0 0.075]), bwperim(max(NucleiMask,[],3)), [1 0 0]);
    PreviewHoechst = imoverlay2(PreviewHoechst, BarMask, [1 1 1]);
    % imtool(PreviewHoechst)
    
    PreviewHoechstHigh = imoverlay2(imadjust(max(ch3,[],3),[0 0.075]), bwperim(max(NucMaskHigh,[],3)), [1 0 0]);
    PreviewHoechstHigh = imoverlay2(PreviewHoechstHigh, BarMask, [1 1 1]);
    % imtool(PreviewHoechstHigh)
    
    PreviewHoechstAlive = imoverlay2(imadjust(max(ch3,[],3),[0 0.075]), bwperim(max(NucMaskAlive,[],3)), [1 0 0]);
    PreviewHoechstAlive = imoverlay2(PreviewHoechstAlive, BarMask, [1 1 1]);
    % imtool(PreviewHoechstAlive)
    
    PreviewTH = imoverlay2(imadjust(max(ch1,[],3),[0 0.05]), bwperim(max(THMask,[],3)), [1 0 0]);
    PreviewTH = imoverlay2(PreviewTH, BarMask, [1 1 1]);
    %imtool(PreviewTH)
  
    PreviewTUJ1 = imoverlay2(imadjust(max(ch4, [], 3), [0 0.09]), bwperim(max(TUJ1Mask,[],3)), [0 0 1]);
    PreviewTUJ1 = imoverlay2(PreviewTUJ1, BarMask, [1 1 1]);
    %imtool(PreviewTUJ1)
 
    PreviewMAP2 = imoverlay2(imadjust(max(ch2, [], 3), [0 0.08]), bwperim(max(MAP2Mask,[],3)), [0 0 1]);
    PreviewMAP2 = imoverlay2(PreviewMAP2, BarMask, [1 1 1]);
    %imtool(PreviewMAP2)
      
    % create empty chanel
    chEmpty = zeros(size(ch1),'uint16');
    

    
    % only Hoechst
    PreviewHoechstAlone = cat (3, imadjust(max(chEmpty,[],3), [0 ; 1],[0 ; 1]), imadjust(max(chEmpty,[],3),[0 ; 1],[0 ; 1]), imadjust(max(ch3,[],3),[0 ; 0.2],[0 ; 1]));
      PreviewHoechstAlone = imoverlay2(PreviewHoechstAlone, BarMask, [1 1 1]);
      %it(PreviewHoechstAlone)
     
    % only TH
    PreviewTHAlone = cat (3, imadjust(max(chEmpty,[],3), [0 ; 1],[0 ; 1]), imadjust(max(ch1,[],3),[0 ; 0.05],[0 ; 1]), imadjust(max(chEmpty,[],3),[0 ; 1],[0 ; 1]));
      PreviewTHAlone = imoverlay2(PreviewTHAlone, BarMask, [1 1 1]);
      %it(PreviewTHAlone)
    
    % only TUJ1
    PreviewTUJ1Alone = cat (3, imadjust(max(ch4,[],3),[0 ; 0.075],[0 ; 1]), imadjust(max(chEmpty,[],3), [0 ; 1],[0 ; 1]), imadjust(max(chEmpty,[],3),[0 ; 1],[0 ; 1]));
       PreviewTUJ1Alone = imoverlay2(PreviewTUJ1Alone, BarMask, [1 1 1]);
       %it(PreviewTUJ1Alone)
    
    % only MAP2
     PreviewMAP2Alone = cat (3, imadjust(max(ch2,[],3), [0 ; 0.09],[0 ; 1]), imadjust(max(chEmpty,[],3),[0 ; 1],[0 ; 1]), imadjust(max(chEmpty,[],3),[0 ; 1],[0 ; 1]));
       PreviewMAP2Alone = imoverlay2(PreviewMAP2Alone, BarMask, [1 1 1]);
       %it(PreviewMAP2Alone)
        
    % combined TH MAP2 TUJ1
    PreviewCombined = cat (3, imadjust(max(ch4,[],3), [0 ; 0.075],[0 ; 1]), imadjust(max(ch1,[],3),[0 ; 0.05],[0 ; 1]), imadjust(max(chEmpty,[],3),[0 ; 0.5],[0 ; 1]));
     PreviewCombined = imoverlay2(PreviewCombined, BarMask, [1 1 1]);
    %it(PreviewCombined)
    
    
    IdentityString = [Label.AreaName{:}, '_Idx_', num2str(Label.Idx)];
    imwrite(PreviewTH, [PreviewPath, filesep, IdentityString, '_', 'TH', '.png'])
    imwrite(PreviewHoechst, [PreviewPath, filesep, IdentityString, '_', 'Hoechst', '.png'])
    imwrite(PreviewHoechstHigh, [PreviewPath, filesep, IdentityString, '_', 'HoechstHigh', '.png'])
    imwrite(PreviewHoechstAlive, [PreviewPath, filesep, IdentityString, '_', 'HoechstAlive', '.png'])
    imwrite(PreviewMAP2, [PreviewPath, filesep, IdentityString, '_', 'MAP2', '.png']) 
    imwrite(PreviewTUJ1, [PreviewPath, filesep, IdentityString, '_', 'TUJ1', '.png'])
    
    imwrite(PreviewHoechstAlone, [PreviewPath, filesep, IdentityString, '_', 'PreviewHoechstAlone', '.png'])
    imwrite(PreviewTHAlone, [PreviewPath, filesep, IdentityString, '_', 'PreviewTHAlone', '.png'])
    imwrite(PreviewTUJ1Alone, [PreviewPath, filesep, IdentityString, '_', 'PreviewTUJ1Alone', '.png'])
    imwrite(PreviewMAP2Alone, [PreviewPath, filesep, IdentityString, '_', 'PreviewMAP2Alone', '.png'])
    imwrite(PreviewCombined, [PreviewPath, filesep, IdentityString, '_', 'PreviewCombined', '.png'])
     
  
    %% Feature extraction
    
    ObjectsThisOrganoid = table();
    ObjectsThisOrganoid.LabelIdx = {Label.Idx};
    ObjectsThisOrganoid.AreaName = {Label.AreaName};
    
    ObjectsThisOrganoid.NucMaskSum = sum(NucleiMask(:));
    ObjectsThisOrganoid.NucMaskHigh = sum(NucMaskHigh(:));
    ObjectsThisOrganoid.NucMaskAlive = sum(NucMaskAlive(:));
    
    ObjectsThisOrganoid.TUJ1MaskSum = sum(TUJ1Mask(:));
    ObjectsThisOrganoid.TUJ1ByNuc = sum(TUJ1Mask(:)) / sum(NucleiMask(:)); 
    ObjectsThisOrganoid.TUJ1ByNucAlive = sum(TUJ1Mask(:)) / sum(NucMaskAlive(:));
    
    ObjectsThisOrganoid.THMaskSum = sum(THMask(:));
    ObjectsThisOrganoid.THByNuc = sum(THMask(:)) / sum(NucleiMask(:));
    ObjectsThisOrganoid.THByNucAlive = sum(THMask(:)) / sum(NucMaskAlive(:));
    
    ObjectsThisOrganoid.MAP2MaskSum = sum(MAP2Mask(:));
    ObjectsThisOrganoid.MAP2ByNuc = sum(MAP2Mask(:)) / sum(NucleiMask(:));
    ObjectsThisOrganoid.MAP2ByNucAlive = sum(MAP2Mask(:)) / sum(NucMaskAlive(:));
    
% %     ObjectsThisOrganoid.THFragmentation = sum(SurfaceTH(:)) / sum(THMask(:));
% %     ObjectsThisOrganoid.SkelTH = sum(skelTH(:));
% %     ObjectsThisOrganoid.Nodes = size(nodeTH, 2);
% %     ObjectsThisOrganoid.Links = size(linkTH, 2);
% %     ObjectsThisOrganoid.THPercent = THPercent;
% %    
    % combinations 
    ObjectsThisOrganoid.MaskTHTUJ1 = sum(MaskTHTUJ1(:));
    ObjectsThisOrganoid.MaskTUJ1MAP2 = sum(MaskTUJ1MAP2(:));
    ObjectsThisOrganoid.MaskMAP2TH = sum(MaskMAP2TH(:));
    
      
end

