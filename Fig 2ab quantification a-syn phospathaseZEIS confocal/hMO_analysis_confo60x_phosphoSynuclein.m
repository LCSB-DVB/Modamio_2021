%% Prepare Matlab


% remove paths added during startup that create problems 
rmpath(genpath('S:\Libraries\hcsforge\CV8000'))
disp('Removed CV8000')
rmpath(genpath('S:\HCS_Platform\Scripts_Repository\Retro'))
disp('Removed Retro from S:\HCS_Platform\Scripts_Repository\Retro')


% code adapted from SB and JJ 
%% Document script

% % 30d
% SavePathMainDirectory = 'S:\HCS_Platform\Data\JenniferModamio\2021\JM20210327_confocal_synPhospho\20210326_30d_phosphatase_synT_synP_Map2\';
% % 70d
 SavePathMainDirectory = '\\atlas\LCSB_HCS\HCS_Platform\Data\JenniferModamio\2021\JM20210327_confocal_synPhospho\20210402_70d_phosphatase_synT_synP_Map2\';


AnalysisTimeStamp = datestr(now, 'yyyymmdd_HHMMSS');
SavePath = [SavePathMainDirectory, 'Analysis_', AnalysisTimeStamp];
mkdir(SavePath)
FileNameShort = mfilename;
%newbackup = sprintf('%s_log.m',[SavePath, '\', FileNameShort]);
%FileNameAndLocation = mfilename('fullpath');
%currentfile = strcat(FileNameAndLocation, '.m');
%copyfile(currentfile,newbackup);
PreviewPath = [SavePath, filesep, 'Previews'];
mkdir(PreviewPath);
Version = version();
save([SavePath filesep 'MatlabVersion.mat'], 'Version')
    
    
   
%% Document Matlab state
%Version = ver;
%MegatronPath = pwd;

%% Analysis
%Objects = {};
ObjectsAll = [];

% % % 30d
%files = dirrec('S:\HCS_Platform\Data\JenniferModamio\2021\JM20210327_confocal_synPhospho\20210326_30d_phosphatase_synT_synP_Map2\', '.czi')';
% % % 70d
 files = dirrec('\\atlas\LCSB_HCS\HCS_Platform\Data\JenniferModamio\2021\JM20210327_confocal_synPhospho\20210402_70d_phosphatase_synT_synP_Map2\', '.czi')';



 Name = regexp(files, '.*\\(.*)', 'tokens'); 
% Name = Name{:};

for s = 1:size(files,1)

    
    sample = files{s};
    data = bfopen(sample); 

    %% Extract metadata from omeXML
    omeMeta = data{1, 4};
    SizeX = omeMeta.getPixelsSizeX(0).getValue(); % image width, pixels
    SizeY = omeMeta.getPixelsSizeY(0).getValue(); % image height, pixels
    SizeZ = omeMeta.getPixelsSizeZ(0).getValue(); % number of Z slices

     SizeC = omeMeta.getPixelsSizeC(0).getValue(); % number of channels
    voxelSizeX = double(omeMeta.getPixelsPhysicalSizeX(0).value(ome.units.UNITS.MICROM)); % in µm
    voxelSizeY = double(omeMeta.getPixelsPhysicalSizeY(0).value(ome.units.UNITS.MICROM)); % in µm


    Fluors = {};
    for chID = 1:SizeC
        Fluors{chID} = char(omeMeta.getChannelFluor(0,chID-1))
    end
    Fluors = {'MAP2', 'SynPhospo', 'Hoechst', 'SynTotal'};

    % Load volume data
    ch1 = cat(3, data{1,1}{(1:SizeC:size(data{1,1},1)),1}); % MAP2  % vol(ch1, 0, 1000)
    ch2 = cat(3, data{1,1}{(2:SizeC:size(data{1,1},1)),1}); % SynTotal    % vol(ch2, 0, 2500)
    ch3 = cat(3, data{1,1}{(3:SizeC:size(data{1,1},1)),1}); %    SynPhospo     vol(ch3, 0, 1000)
    ch4 = cat(3, data{1,1}{(4:SizeC:size(data{1,1},1)),1}); % Hoechst   % vol(ch4, 0, 100)
    
     %tic
    psf_Red = load('psf_red.mat');
    psf_Red = psf_Red.psf; %vol(psf_Red)
    [RedDeconvolvedIm, PSF] = deconvblind(ch2, psf_Red, 10);
    %toc 5 minutes
    
    psf_Green = load('psf_green.mat');
    psf_Green = psf_Green.psf; %vol(GreenDeconvolvedIm)
    [GreenDeconvolvedIm, PSF] = deconvblind(ch3, psf_Green, 10); % Matlab internal function

    %% Initialize variables
    NucleiMask = [];
    SynPhospoMask = [];
    SynTotalMask = [];
    MAP2Mask = [];
    
    %% Segment nuclei (ch4)   
     
    ch4BlurSmall = imfilter(ch4, fspecial('gaussian', 10, 1), 'symmetric');%vol(ch4BlurSmall)
    ch4BlurBig = imfilter(ch4, fspecial('gaussian', 300, 100), 'symmetric');%vol(ch3BlurBig) %%kind of flatfield corraction, to account for different bk in the pic
    ch4DoG = ch4BlurSmall - ch4BlurBig; % vol(ch4DoG, 0, 30)
    
    NucleiMaskDoG = ch4DoG > 300; %vol(NucleiMaskDoG)
    NucleiMask = bwareaopen(NucleiMaskDoG, 20); %vol(NucleiMask)

    %% MAP2 (ch1) vol(ch1)
    ch1MedFilt = [];
    for p = 1:SizeZ
        ch1MedFilt(:,:,p) = medfilt2(ch1(:,:,p));
    end
    ch1MedFilt(:,:,p) = medfilt2(ch1(:,:,p)); %vol(ch1MedFilt)
    MAP2Mask = ch1MedFilt > 500; %vol(MAP2Mask) 
    MAP2Mask = bwareaopen(MAP2Mask, 5);
   
  
    %% SynTotal (Presynaptic structures)
      
    SynTotalMed = medfilt3(RedDeconvolvedIm); %vol(SynTotalMed, 0, 100, 'hot') %vol(SynTotalMed)  
    SynTotalMed_Mask = SynTotalMed > 300; %vol(SynTotalMed_Mask)
    SynTotalMask = SynTotalMed_Mask & MAP2Mask;
    %SynTotalMask = f_RemoveBigObjects(SynTotalMask,10000); 
    %vol(SynTotalMask)
    
    %% Synuclein filament
    SynPhospoMed = medfilt3(GreenDeconvolvedIm); %vol(SynPhospoMed, 0, 30, 'hot')
    SynPhospoMed_Mask = SynPhospoMed > 300; %vol(SynPhospoMed_Mask) 
    SynPhospoMask = SynPhospoMed_Mask & MAP2Mask;
    %SynPhospoMask = f_RemoveBigObjects(SynPhospoMask,10000); 
    %vol(SynPhospoMask, 0, 1)
   
             
    %% SynTotal co-localization 
    
    SynPhospoSynTotalMap2Mask = SynPhospoMask & SynTotalMask ;  %vol(SynPhospoSynTotalMap2Mask)
    SynPhospoNoNSynTotalMap2Mask = SynPhospoMask & ~SynTotalMask ; %vol(SynPhospoNoNSynTotalMap2Mask)
    SynTotalNoNSynPhospoMap2Mask = SynTotalMask & ~SynPhospoMask ; %vol(SynTotalNoNSynPhospoMap2Mask)
     
    %% Previews 
      % Scalebar
      
%     imSize = size(ch1);
%     [BarMask, BarCenter] = f_barMask(20, voxelSizeX, imSize, imSize(1)-100, 100, 10);
     
    % Scalebar
    imSize = [size(ch1, 1), size(ch1, 2)];
    [BarMask, BarCenter] = f_barMask(200, 0.42, imSize, imSize(1)-200, 200, 25);
    
    ImHoechst = imadjust(max(ch4,[],3), [0,0.09], [0,1]); %it(ImHoechst)
    ImSynTotal = imadjust(max(ch2,[],3), [0,0.09], [0,1]); %it(ImSynTotal)
    ImSynPhospo = imadjust(max(ch3,[],3), [0,0.09], [0,1]); %it(ImSynPhospo)
    ImMAP2 = imadjust(max(ch1,[],3), [0,0.5], [0,1]); %it(ImMAP2)
    
    imwrite(ImHoechst, [PreviewPath, filesep, Name{s,1}{1,1}{1,1},'ImHoechst.png']);
    imwrite(ImSynTotal, [PreviewPath, filesep, Name{s,1}{1,1}{1,1},'ImSynTotal.png']);
    imwrite(ImSynPhospo, [PreviewPath, filesep, Name{s,1}{1,1}{1,1},'ImSynPhospo.png']);
    imwrite(ImMAP2, [PreviewPath,filesep, Name{s,1}{1,1}{1,1},'ImMAP2.png']);
     
     
    PreviewNuclei= imoverlay2(imadjust(max(ch4,[],3), [0,0.04], [0,1]), bwperim(max(NucleiMask,[],3)), [1 0 0]);
    PreviewNuclei = imoverlay2(PreviewNuclei, BarMask, [1 1 1]);
    % it(PreviewNuclei)
    PreviewSynTotal= imoverlay2(imadjust(max(ch2,[],3), [0,0.04], [0,1]), bwperim(max(SynTotalMask,[],3)), [1 0 0]);
    PreviewSynTotal = imoverlay2(PreviewSynTotal, BarMask, [1 1 1]);
    % it(PreviewSynTotal)
    PreviewSynPhospo= imoverlay2(imadjust(max(ch3,[],3), [0,0.04], [0,1]), bwperim(max(SynPhospoMask,[],3)), [1 0 0]);
    PreviewSynPhospo = imoverlay2(PreviewSynPhospo, BarMask, [1 1 1]);
    % it(PreviewSynPhospo)
    PreviewMAP2 =imoverlay2(imadjust(max(ch1,[],3), [0,0.05], [0,1]), bwperim(max(MAP2Mask,[],3)), [1 0 0]);
    PreviewMAP2 = imoverlay2(PreviewMAP2, BarMask, [1 1 1]);
    % it(PreviewMAP2)

    PreviewAll= imoverlay2(imadjust(max(ch1,[],3), [0,0.02], [0,1]), bwperim(max(MAP2Mask,[],3)), [0 1 1]);%blue
    PreviewAll = imoverlay2(PreviewAll, bwperim(max(SynTotalMask,[],3)),  [1 0 0]); %red 
    PreviewAll = imoverlay2(PreviewAll, bwperim(max(SynPhospoMask,[],3)),[0 1 0]); %green 
    % it(PreviewAll)
    
    PreviewSyn= imoverlay2(imadjust(max(ch2,[],3), [0,0.04], [0,1]), bwperim(max(SynTotalMask,[],3)), [1 0 0]);%red
    PreviewSyn = imoverlay2(PreviewSyn, bwperim(max(SynPhospoMask,[],3)),[0 1 0]); %green;
    % it(PreviewSyn)
      
    imwrite(PreviewNuclei, [PreviewPath, filesep, Name{s,1}{1,1}{1,1},'_PreviewNuclei.png']); 
    imwrite(PreviewSynTotal, [PreviewPath, filesep, Name{s,1}{1,1}{1,1},'_PreviewSynTotal.png']);
    imwrite(PreviewSynPhospo, [PreviewPath, filesep, Name{s,1}{1,1}{1,1},'_PreviewSynPhospo95.png']);
    imwrite(PreviewMAP2, [PreviewPath, filesep, Name{s,1}{1,1}{1,1},'_PreviewMAP2.png']);
    imwrite(PreviewAll, [PreviewPath,filesep, Name{s,1}{1,1}{1,1},'_PreviewAll.png']);
    imwrite(PreviewSyn, [PreviewPath,filesep, Name{s,1}{1,1}{1,1},'_PreviewSyn.png']);

    
    %% Collect data
    ObjectsThisField = table();
    ObjectsThisField.sample = {sample};
    ObjectsThisField.NeuroVolume = sum(MAP2Mask(:));
    ObjectsThisField.NucVolume = sum(NucleiMask(:));
    ObjectsThisField.SynTotalVolume = sum(SynTotalMask(:));
    ObjectsThisField.SynPhosphoVolume = sum(SynPhospoMask(:));
    
    ObjectsThisField.SynPhospoSynTotalMap2Mask = sum(SynPhospoSynTotalMap2Mask(:));
    ObjectsThisField.SynPhospoNoNSynTotalMap2Mask = sum(SynPhospoNoNSynTotalMap2Mask(:));
    ObjectsThisField.SynTotalNoNSynPhospoMap2Mask = sum(SynTotalNoNSynPhospoMap2Mask(:));
    
    if s == 1
        ObjectsAll = ObjectsThisField;
    else
        ObjectsAll = [ObjectsAll; ObjectsThisField];
    end
   
end

save([SavePath, filesep, 'workspace.mat'], 'ObjectsAll', 'Version')
writetable(ObjectsAll, [SavePath, filesep, 'Objects.xls'])
    
    
    
    
    
    
    
    
    
    
    
    