%% Prepare Matlab
clear
clc

% remove paths creating conflict
% remove paths added during startup that create problems 
rmpath(genpath('S:\Libraries\hcsforge\CV8000'))
disp('Removed CV8000')
rmpath(genpath('S:\HCS_Platform\Scripts_Repository\Retro'))
disp('Removed Retro from S:\HCS_Platform\Scripts_Repository\Retro')

% code adapted from SB and JJ 
%% Document script


%SavePathMainDirectory = 'S:\HCS_Platform\Data\JenniferModamio\2020\hMO_characterisation\JM20200714_SNCA_Synapse_PDS95_Synapthophysin_MAP2\synapse_15d\';
%SavePathMainDirectory = 'S:\HCS_Platform\Data\JenniferModamio\2020\hMO_characterisation\JM20200714_SNCA_Synapse_PDS95_Synapthophysin_MAP2\synapse_30d\';
SavePathMainDirectory = 'S:\HCS_Platform\Data\JenniferModamio\2020\hMO_characterisation\JM20200714_SNCA_Synapse_PDS95_Synapthophysin_MAP2\synapse_70d\';


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
%files = dirrec('X:\groups\schwamborn\jennifer.modamio\2020\20200621_synapse\30day diff hMO\', '.czi')';
% files = dirrec('X:\groups\schwamborn\jennifer.modamio\2020\20200621_synapse\70day diff hMO\', '.czi')';
% files = dirrec('X:\groups\schwamborn\jennifer.modamio\2020\20200621_synapse\15day diff hMO\', '.czi')';

% analysis from hcs2 
%files = dirrec('S:\HCS_Platform\Data\JenniferModamio\2020\Images\synapse images analysis 2020\15day diff hMO\', '.czi')';
%files = dirrec('S:\HCS_Platform\Data\JenniferModamio\2020\Images\synapse images analysis 2020\30day diff hMO\', '.czi')';
files = dirrec('S:\HCS_Platform\Data\JenniferModamio\2020\Images\synapse images analysis 2020\70day diff hMO\', '.czi')';


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
%     SizeT = omeMeta.getPixelsSizeT(0).getValue(); % number of timepoints
     SizeC = omeMeta.getPixelsSizeC(0).getValue(); % number of channels
    voxelSizeX = double(omeMeta.getPixelsPhysicalSizeX(0).value(ome.units.UNITS.MICROM)); % in µm
    voxelSizeY = double(omeMeta.getPixelsPhysicalSizeY(0).value(ome.units.UNITS.MICROM)); % in µm
%     voxelSizeZ = double(omeMeta.getPixelsPhysicalSizeZ(0).value(ome.units.UNITS.MICROM)); % in µm

    Fluors = {};
    for chID = 1:SizeC
        Fluors{chID} = char(omeMeta.getChannelFluor(0,chID-1))
    end
    Fluors = {'MAP2', 'PSD95', 'Hoechst', 'Synaptophysin'};

    % Load volume data
    ch1 = cat(3, data{1,1}{(1:SizeC:size(data{1,1},1)),1}); % MAP2  % vol(ch1, 0, 1000)
    ch2 = cat(3, data{1,1}{(2:SizeC:size(data{1,1},1)),1}); % Synaptophysin    % vol(ch2, 0, 2500)
    ch3 = cat(3, data{1,1}{(3:SizeC:size(data{1,1},1)),1}); %    PSD95     vol(ch3, 0, 1000)
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
    PSD95Mask = [];
    SynaptophysinMask = [];
    MAP2Mask = [];
    
    %% Segment nuclei (ch4)   
    ch4BlurSmall = imfilter(ch4, fspecial('gaussian', 10, 1), 'symmetric');%vol(ch4BlurSmall)
    ch4BlurBig = imfilter(ch4, fspecial('gaussian', 300, 100), 'symmetric');%vol(ch3BlurBig) %%kind of flatfield corraction, to account for different bk in the pic
    ch4DoG = ch4BlurSmall - ch4BlurBig; % vol(ch4DoG, 0, 30)
    
    NucleiMaskDoG = ch4DoG>5; %vol(NucleiMaskDoG)
    NucSplitter = imfilter(ch4, fspecial('log', 10, 2), 'symmetric'); %vol(NucSplitter, 0, 5)
    % figure; surf(fspecial('log', 10, 2))
    NucSplitter = NucSplitter > 0; %catch the pixel for identifying nuclear splitter (this is the valley of the previous figure)
    %vol(NucSplitter, 0, 1)
    NucSplitter = bwareaopen(NucSplitter, 500); %get ride of nucleoli for example
    NucleiMask = NucleiMaskDoG & ~NucSplitter; %vol(NucleiMask)
    
    %% MAP2 (ch1) vol(ch1)
    ch1MedFilt = [];
    for p = 1:SizeZ
        ch1MedFilt(:,:,p) = medfilt2(ch1(:,:,p));
    end
    ch1MedFilt(:,:,p) = medfilt2(ch1(:,:,p)); %vol(ch1MedFilt)
    MAP2Mask = ch1MedFilt > 500; %vol(MAP2Mask) 
    MAP2Mask = bwareaopen(MAP2Mask, 20);
   
    MAP2MaskPerim = imdilate(MAP2Mask, strel('disk',4));  %vol(MAP2MaskPerim)
    %% Synaptophysin (Presynaptic structures)
      
    %vol(GreenDeconvolvedIm, 0, 10000, 'hot')
    %vol(ch2, 0, 1000, 'hot')
    SynaptophysinMed = medfilt3(RedDeconvolvedIm); %vol(SynaptophysinMed, 0, 100, 'hot') %vol(SynaptophysinMed)  
    SynaptophysinMed_Mask = SynaptophysinMed > 400; %vol(SynaptophysinMed_Mask)
    SynaptophysinDoG = imfilter(RedDeconvolvedIm, fspecial('gaussian', 11, 1), 'symmetric') - imfilter(RedDeconvolvedIm, fspecial('gaussian', 11, 3), 'symmetric');
    %vol(SynaptophysinDoG, 0,50,'hot')vol(SynaptophysinMask)
    SynaptophysinMask = SynaptophysinDoG > 100;%
    SynaptophysinMask = SynaptophysinMask & SynaptophysinMed_Mask;

    SynaptophysinMask = SynaptophysinMask & MAP2Mask & ~NucleiMask;
    %SynaptophysinMask = bwareaopen(SynaptophysinMask, 27); % 27 neighbouhood idea
    SynaptophysinMask = f_RemoveBigObjects(SynaptophysinMask, 2000);
    %vol(SynaptophysinMask)
    
    %% PSD95 - post-synaptic compartment
    PSDMed = medfilt3(GreenDeconvolvedIm); %vol(PSDMed, 0, 30, 'hot')
    PSDMed_Mask = PSDMed > 500; %vol(PSDMed_Mask)
    PSDDoG = imfilter(GreenDeconvolvedIm, fspecial('gaussian', 11, 1), 'symmetric') - imfilter(GreenDeconvolvedIm, fspecial('gaussian', 11, 3), 'symmetric');
    %vol(PSDDoG, 0,10,'hot')
    PSDMask = PSDDoG > 150; %vol(PSDMask)
    PSDMask = PSDMask & MAP2Mask & ~NucleiMask; %vol(PSDMask, 0, 1)
    %PSDMask = bwareaopen(PSDMask,27 ); % 27 neighbouhood idea
    PSDMask = f_RemoveBigObjects(PSDMask,2000); %vol(PSDMask, 0, 1)
    
%     
%     PSDMaskMap2Perim = PSDMask & MAP2MaskPerim & ~NucleiMask;
%     PSDMaskMap2Perim = f_RemoveBigObjects(MAP2MaskPerim, 1000); %vol(PSDMaskMap2Perim, 0, 1)
%     
    %% Synaptophysin neighbourhood
    SynNeighbourhood = imdilate(SynaptophysinMask, strel('disk',4));
    % vol(SynaptophysinMask + 2* SynNeighbourhood, 0,3,'jet')
    PSDNeighbourhood = imreconstruct(SynNeighbourhood, PSDMask);
    % vol(SynaptophysinMask + 2*PSDNeighbourhood,0,3,'jet')
    SynInSynapses = imreconstruct(imdilate(PSDNeighbourhood, strel('disk',5)), SynaptophysinMask);
    % vol(SynInSynapses + 2*PSDNeighbourhood,0,3,'jet')
               
    %% Previews 
      % Scalebar
%     imSize = size(ch1);
%     [BarMask, BarCenter] = f_barMask(20, voxelSizeX, imSize, imSize(1)-100, 100, 10);
%     
      % Scalebar
    imSize = [size(ch1, 1), size(ch1, 2)];
    [BarMask, BarCenter] = f_barMask(200, 0.42, imSize, imSize(1)-200, 200, 25);
    
    ImHoechst = imadjust(max(ch4,[],3), [0,0.09], [0,1]); %it(ImHoechst)
    ImSynaptophysin = imadjust(max(ch2,[],3), [0,0.4], [0,1]); %it(ImSynaptophysin)
    ImPSD = imadjust(max(ch3,[],3), [0,0.3], [0,1]); %it(ImPSD)
    ImMAP2 = imadjust(max(ch1,[],3), [0,0.5], [0,1]); %it(ImMAP2)
    
    
    PreviewSynaptophysin= imoverlay2(imadjust(max(ch2,[],3), [0,0.02], [0,1]), bwperim(max(SynaptophysinMask,[],3)), [1 0 0]);
    PreviewSynaptophysin = imoverlay2(PreviewSynaptophysin, BarMask, [1 1 1]);
    % it(PreviewSynaptophysin)
    PreviewPSD= imoverlay2(imadjust(max(ch3,[],3), [0,0.02], [0,1]), bwperim(max(PSDMask,[],3)), [1 0 0]);
    PreviewPSD = imoverlay2(PreviewPSD, BarMask, [1 1 1]);
    % it(PreviewPSD)
    PreviewMAP2 =imoverlay2(imadjust(max(ch1,[],3), [0,0.02], [0,1]), bwperim(max(MAP2Mask,[],3)), [1 0 0]);
    PreviewMAP2 = imoverlay2(PreviewMAP2, BarMask, [1 1 1]);
    % it(PreviewMAP2)

    PreviewSynptPSDMAP2= imoverlay2(imadjust(max(ch1,[],3), [0,0.02], [0,1]), bwperim(max(PSDMask,[],3)), [1 0 0]);
    PreviewSynptPSDMAP2 = imoverlay2(PreviewSynptPSDMAP2, bwperim(max(SynaptophysinMask,[],3)), [0 1 0]);
    % it(PreviewMain)
    
    PreviewMAP2Synapses =imoverlay2(imadjust(max(ch1,[],3), [0,0.02], [0,1]), bwperim(max(PSDNeighbourhood,[],3)), [1 0 0]);
    PreviewMAP2Synapses = imoverlay2(PreviewMAP2Synapses, bwperim(max(SynInSynapses,[],3)), [0 1 0]);
   % it(PreviewMAP2Synapses)
   
     imwrite(PreviewSynaptophysin, [PreviewPath, filesep, Name{s,1}{1,1}{1,1},'_PreviewSynaptophysin.png']);
     imwrite(PreviewPSD, [PreviewPath, filesep, Name{s,1}{1,1}{1,1},'_PreviewPSD95.png']);
     imwrite(PreviewMAP2, [PreviewPath, filesep, Name{s,1}{1,1}{1,1},'_PreviewMAP2.png']);
     imwrite(PreviewSynptPSDMAP2, [PreviewPath,filesep, Name{s,1}{1,1}{1,1},'_PreviewSynptPSDMAP2.png']);
     imwrite(PreviewMAP2Synapses, [PreviewPath,filesep, Name{s,1}{1,1}{1,1},'_PreviewMAP2Synapses.png']);
    
   
%      imwrite(PreviewSynaptophysin, [PreviewPath, filesep, fileName{1,s}{1,1}{1,1},'_PreviewSynaptophysin.png']);
%      imwrite(PreviewPSD, [PreviewPath, filesep, fileName{1,s}{1,1}{1,1},'_PreviewPSD95.png']);
%      imwrite(PreviewMAP2, [PreviewPath, filesep, fileName{1,s}{1,1}{1,1},'_PreviewMAP2.png']);
%      imwrite(PreviewSynptPSDMAP2, [PreviewPath,filesep, fileName{1,s}{1,1}{1,1},'_PreviewSynptPSDMAP2.png']);
%      imwrite(PreviewMAP2Synapses, [PreviewPath,filesep, fileName{1,s}{1,1}{1,1},'_PreviewMAP2Synapses.png']);
    

   
    
    %% Collect data
    ObjectsThisField = table();
    ObjectsThisField.sample = {sample};
    ObjectsThisField.SynaptophysinVolume = sum(SynaptophysinMask(:));
    [~, ObjectsThisField.SynaptophysinCount]= bwlabeln(SynaptophysinMask,26);
    ObjectsThisField.PSDVolume = sum(PSDMask(:));
    [~, ObjectsThisField.PSDCount] = bwlabeln(PSDMask, 26);
     ObjectsThisField.SynaptophysinSynapseVolume = sum(SynInSynapses(:));
    [~, ObjectsThisField.SynaptophysinSynapseCount] = bwlabeln(SynInSynapses, 26);
     ObjectsThisField.PSDSynapseVolume = sum(PSDNeighbourhood(:));
    [~, ObjectsThisField.PSDSynapseCount] = bwlabeln(PSDNeighbourhood, 26);
    ObjectsThisField.NeuroVolume = sum(MAP2Mask(:));
    ObjectsThisField.NucVolume = sum(NucleiMask(:));

  
    if s == 1
        ObjectsAll = ObjectsThisField;
    else
        ObjectsAll = [ObjectsAll; ObjectsThisField];
    end
   
end

save([SavePath, filesep, 'workspace.mat'], 'ObjectsAll', 'Version')
writetable(ObjectsAll, [SavePath, filesep, 'Objects.xls'])
    
    
    
    
    
    
    
    
    
    
    
    