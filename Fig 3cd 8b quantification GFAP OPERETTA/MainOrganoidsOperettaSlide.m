%% User inputs
addpath(genpath('S:\Libraries\hcsforge\LCSBMatlabLibrary'))
disp('Loaded LCSB Matlab Library')

clear
clc

SetupMode = 0; % 1 for creating numeric organoid labels OR 0 for linking the final analysis to human labels

%% all samples done % 

%%%%%%% text file defining the samples %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Astrocytes 15 days 
% done SlideLayout = 'JM20191126_DB322_DB317_GFAP_647_S100B_568_TUJ1_488_15dDiff.txt';
% % Astrocytes 30 days 
% SlideLayout = 'JM20191126_DB322_DB317_GFAP_647_S100B_568_TUJ1_488_30dDiff.txt';
% % Astrocytes 50 days 
% SlideLayout = 'JM20191126_DB322_DB317_GFAP_647_S100B_568_TUJ1_488_50dDiff.txt';
% % Astrocytes 70 days 
% done SlideLayout = 'JM20191126_DB322_DB317_GFAP_647_S100B_568_TUJ1_488_70dDiff.txt';
% % Astrocytes 90 days 
% SlideLayout = 'JM20191126_DB322_DB317_GFAP_647_S100B_568_TUJ1_488_90dDiff.txt';
% % Astrocytes 110 days 
% SlideLayout = 'JM20191126_DB322_DB317_GFAP_647_S100B_568_TUJ1_488_110dDiff.txt';


% save path %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % Astrocytes 15 days
 SavePath = 'S:\HCS_Platform\Data\JenniferModamio\2019\JM_2019DB317DB322_analysisDB322DB317_firstDerivation\Astrocytes\JM20191126_DB322_DB317_GFAP_647_S100B_568_TUJ1_488_15dDiff';
% % Astrocytes 30 days
%  done SavePath = 'S:\HCS_Platform\Data\JenniferModamio\2019\JM_2019DB317DB322_analysisDB322DB317_firstDerivation\Astrocytes\JM20191126_DB322_DB317_GFAP_647_S100B_568_TUJ1_488_30dDiff';
% % Astrocytes 50 days
% done SavePath = 'S:\HCS_Platform\Data\JenniferModamio\2019\JM_2019DB317DB322_analysisDB322DB317_firstDerivation\Astrocytes\JM20191126_DB322_DB317_GFAP_647_S100B_568_TUJ1_488_50dDiff';
% % Astrocytes 70 days
% done  SavePath = 'S:\HCS_Platform\Data\JenniferModamio\2019\JM_2019DB317DB322_analysisDB322DB317_firstDerivation\Astrocytes\JM20191126_DB322_DB317_GFAP_647_S100B_568_TUJ1_488_70dDiff';
% % Astrocytes 90 days
%SavePath = 'S:\HCS_Platform\Data\JenniferModamio\2019\JM_2019DB317DB322_analysisDB322DB317_firstDerivation\Astrocytes\JM20191126_DB322_DB317_GFAP_647_S100B_568_TUJ1_488_90dDiff';
% % Astrocytes 110 days
% SavePath = 'S:\HCS_Platform\Data\JenniferModamio\2019\JM_2019DB317DB322_analysisDB322DB317_firstDerivation\Astrocytes\JM20191126_DB322_DB317_GFAP_647_S100B_568_TUJ1_488_110dDiff';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


mkdir(SavePath)
PreviewSavePath = [SavePath, filesep, 'Previews'];
mkdir(PreviewSavePath)

%% Parallel pool control
delete(gcp('nocreate'))
myCluster = parcluster;
Workers = myCluster.NumWorkers;
% parpool(28) %for HURRICANE
parpool(Workers) % for MEGATRON

%% Run mode control
if SetupMode
    RunMode = 0;
else
    RunMode = 1;
end

%% Common part

channelID = 1; % channel to use for overview

% Input information %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % Astrocytes 15 days
InfoTable = readtable('S:\Operetta\OperettaDB_LCSB\JM20191126_DB322_DB317_GFAP_647_S100B_568_TUJ1_488_15dDiff\2aaf1e86-4054-496c-b92c-e0f7c4631950\metadata.csv');
% % Astrocytes 30 days
% InfoTable = readtable('S:\Operetta\OperettaDB_LCSB\JM20191126_DB322_DB317_GFAP_647_S100B_568_TUJ1_488_30dDiff\f6c9fcfe-6e86-4719-b283-06d524dd88db\metadata.csv');
% % Astrocytes 50 days
% InfoTable = readtable('S:\Operetta\OperettaDB_LCSB\JM20191126_DB322_DB317_GFAP_647_S100B_568_TUJ1_488_50dDiff\240ce419-3c80-4b97-b04b-f4b94693321e\metadata.csv');
% % Astrocytes 70 days
% done  InfoTable = readtable('S:\Operetta\OperettaDB_LCSB\JM20191126_DB322_DB317_GFAP_647_S100B_568_TUJ1_488_70dDiff\27a944d3-7812-46f5-85e9-ba75bb034683\metadata.csv');
% % Astrocytes 90 days
% InfoTable = readtable('S:\Operetta\OperettaDB_LCSB\JM20191126_DB322_DB317_GFAP_647_S100B_568_TUJ1_488_90dDiff\d81c9aa2-7fa9-4e59-acfb-ac0ab0d13a58\metadata.csv');
% % Astrocytes 110 days
% InfoTable = readtable('S:\Operetta\OperettaDB_LCSB\JM20191126_DB322_DB317_GFAP_647_S100B_568_TUJ1_488_110dDiff\7742d864-ebd0-437b-8e55-79d892dc3716\metadata.csv');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ChannelNames = unique(InfoTable.Channel);
    Channels = length(unique(InfoTable.Channel));
    Planes = unique(InfoTable.Plane)';
    Timepoints = unique(InfoTable.Timepoint)' + 1;
    [GroupsTable, GroupsIm5DCellArray] = FindGroups(InfoTable); % it(GroupsTable)
    
%% Setup mode

if SetupMode == 1
    
    Preview = CreateLabelHelpPreview(GroupsTable, PreviewSavePath);
    imwrite(Preview, [PreviewSavePath, filesep, 'layout.png'])
    Message = ['The plate Layout has been saved at ', [PreviewSavePath, filesep, 'layout.png'], '. Please save a text file without header, using tab separation, where the first column is the index number as shown in the preview and the second is the area name. Save the text file as SlideLayout_Date.txt in your working directory and set the variable SetupMode to 0 >>> Run'];
    h = msgbox(Message);
    
else    
%% Analysis mode
    
    % Load annotations
    Layout = readtable(SlideLayout)
    Layout.Properties.VariableNames = {'Idx', 'AreaName'};
    
    % Load images and organize in an XYC array
    Groups = unique(GroupsTable(GroupsTable > 0))';
    GroupPreviews = {};
    ObjectsAll = {};
    

    
    for g = Groups % Number of organoids
         XYMosaicCells = {};
        GroupZone = GroupsTable == g;
        [GroupIdxRowVec, GroupIdxColVec] = find(GroupZone); % linear indexes
        Elements = sum(GroupZone(:));
        InfoTablesThisGroup = {};
        for e = 1:Elements % Fields of a given organoid
            for c = 1:Channels
            InfoTableThisField = GroupsIm5DCellArray{GroupIdxRowVec(e), GroupIdxColVec(e)};
            InfoTablesThisGroup{e} = InfoTableThisField;
            InfoTableThisChannel = InfoTableThisField(strcmp(InfoTableThisField.Channel, ChannelNames{c}), :);
                clear Im4D
                for t = Timepoints
                    for p = Planes
                        InfoTableThisChannelThisPlane = InfoTableThisChannel(InfoTableThisChannel.Plane == p, :);
                        ImPathThisPlane = InfoTableThisChannelThisPlane.Path{:};   
                        Im4D(:,:,t,p) = imread(ImPathThisPlane); % it(Im4D(:,:,t,p))
                    end
                end
               XYMosaicCells{c}{GroupIdxRowVec(e), GroupIdxColVec(e)} = Im4D; % Operetta counterpart of XYmosaicCells for Opera
            end
        end

        InfoTableThisGroup = vertcat(InfoTablesThisGroup{:});

        %% Remove empty cells
        XYMosaicCells = cellfun(@(x) GroupClipper(x),  XYMosaicCells, 'UniformOutput', false);

        %% Stitch
        XYmosaicContourCell = cellfun(@(x) stdfilt(x, ones(3)), XYMosaicCells{1}, 'UniformOutput', false);
        XPositions = unique(InfoTableThisGroup.PositionX); % m
        YPositions = unique(InfoTableThisGroup.PositionY); % m
        ResolutionXY = 675 / 1360; % um per pixel
        MaxPixelDrift = 30;
        PreviewChannel = 1;
        ShowProgress = 0;
        [CroppedMosaic, StitchedIm] = f_stitching_operetta(XYMosaicCells, XYmosaicContourCell, XPositions, YPositions, ResolutionXY, MaxPixelDrift, PreviewChannel, ShowProgress);
        GroupPreviews{g} = max(CroppedMosaic{channelID},[],3); %it(GroupPreviews{g})
        
        %% Image analysis
        Label = Layout(g,:);
        try
            ObjectsThisOrganoid = f_ImageAnalysisPerOperettaOrganoid_cell_count(Label, CroppedMosaic{1}, CroppedMosaic{2}, CroppedMosaic{3}, CroppedMosaic{4}, ChannelNames, PreviewSavePath);
        catch
            Errors{g} = 'Image analysis failed';
            continue % next group g
        end
            ObjectsAll{g} = ObjectsThisOrganoid;
    end
    
    Objects = vertcat(ObjectsAll{:});
    save([SavePath, filesep, 'Objects.mat'], 'Objects');
    writetable(Objects, [SavePath, filesep, 'Objects.csv'])
    writetable(Objects, [SavePath, filesep, 'Objects.xlsx'])

    %% Preview of the whole slide

    SizeSingleIm = size(XYMosaicCells{1}{1,1});
    SizeSingleIm = SizeSingleIm(1:2);
    RescanGridSize = size(GroupsTable);
    GreatPreview = zeros(SizeSingleIm(1)*RescanGridSize(1), SizeSingleIm(2)*RescanGridSize(2), 'uint16');
    ImHeight = SizeSingleIm(1);
    ImWidth = SizeSingleIm(2);
    StartRCell = {};
    StartCCell = {};

    for g = Groups
        g
        StitchedGroupSize = size(GroupPreviews{g});
        ZoneNow = GroupsTable == g;
        [R,C] = find(ZoneNow)
        StartR = min(R);
        StartC = min(C);
        StartRPixel = ((StartR-1) * ImHeight) + 1;
        %EndRPixel = StartRPixel + (3 * ImHeight) - 1;
        EndRPixel = StartRPixel + StitchedGroupSize(1) - 1;
        StartCPixel = ((StartC-1) * ImWidth) + 1;
        %EndCPixel = StartCPixel + (3 * ImWidth) - 1;
        EndCPixel = StartCPixel + StitchedGroupSize(2) - 1;
        GreatPreview(StartRPixel:EndRPixel, StartCPixel:EndCPixel) = GroupPreviews{g};
        StartRCell{g} = StartRPixel;
        StartCCell{g} = StartCPixel;
    end

    Zoomfactor = 50;
    %GreatPreviewResized = imresize(imadjust(GreatPreview), 1/Zoomfactor);
    GreatPreviewResized = imresize(imadjust(GreatPreview, [0 0.02], [0 1]), 1/Zoomfactor);

    for g = Groups
        GreatPreviewResized = insertText(GreatPreviewResized, [round(StartCCell{g}/Zoomfactor), round(StartRCell{g}/Zoomfactor)], num2str(g), 'FontSize', 12, 'BoxColor', 'red', 'TextColor', 'white');
    end
    
    %imtool(GreatPreview)
    %imtool(GreatPreviewResized)
    imwrite(GreatPreviewResized, [SavePath, filesep, 'GreatPreview.png'])
    
    % it(GreatPreviewResized)
    % save([SavePath, filesep, 'WorkspaceIncludingObjects.mat'])
    
end




