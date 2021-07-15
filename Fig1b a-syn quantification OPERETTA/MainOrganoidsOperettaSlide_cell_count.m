%% JM 13.03.2019
%S:\HCS_Platform\Scripts_Repository\JenniferModamio\2019SNCA\JM20190315_SNCAsynuclein\ScriptSyn20190709

%% stainings: 
% 421 total synuclein ms 1:1000 568
% 439 syn filament rb 1:5000 488
% 111 TH chk 1:1000 647 

% the stainig for TH is commonly very difficult to detec and quantify.
% focus only in synuclein (and common with TH), but not TH fragmentation
% and skeleton

%% User inputs
clear
clc

SetupMode = 0; % 1 for creating numeric organoid labels OR 0 for linking the final analysis to human labels

%% first batch of organoids 

% DB322
% SlideLayout = 'JM20190316_DB322_slideNum6_synFil488_syn2568_TH647.txt'; 
% % DB317
% SlideLayout = 'JM20190316_DB317_slideNum5_synFil488_syn2568_TH647.txt'; 
% DB322 slide 12.1
% SlideLayout = 'JM20190402_DB322_slideNum12_SynFilamnt488_SynTotal568_TH647_25_5.txt'; 
% DB322 slide 12.2
% SlideLayout = 'JM20190402_DB322_slideNum12_SynFilamnt488_SynTotal568_TH647_40_6.txt'; 
% % DB317
% SlideLayout = 'JM20190402_DB317_slideNum11_synFil488_syn2568_TH647_26_53um_all.txt'; 

%% batch test old organoids 

% DB322 & 317 70d old (synuclein-second colum each condition)
% SlideLayout = 'JM20190410_DB317_th_syn_stem_322_th_syn_stem_p3_70d.txt'; 

%% batch from 27/05/19 approx 

% DB322 & 317 30d old 
% SlideLayout = 'JM20190529_DB322_DB317_p10p12p14_30dDiff_synFilamnt488_synTotal568_TH647_.txt'; 

% DB322 & 317 50d old 
%SlideLayout = 'JM20190528_DB322_DB317_p8p10p11p12_50dDiff_synFilamnt488_synTotal568_TH647.txt'; 

% DB322 & 317 70d old 
% SlideLayout = 'JM20190528_DB322_DB317_p6p10p11p12_70dDiff_synFilamnt488_synTotal568_TH647.txt'; 

%% batch from 20/07/19 approx 
SlideLayout = 'JM20190727_322_317_p10_p11_p12_90dDiff_filmntsyn_488_TotalSyn_56.txt'; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% savepath %%%%%%%%%%%%%%%%%%

%% first batch of organoids 
% DB322
% SavePath = 'S:\HCS_Platform\Data\JenniferModamio\2019DB317DB322_analysis\Synuclein\ScriptSyn20190709_JM20190316_DB322_slideNum6_synFil488_syn2568_TH647_testHighSynuclein4level'; % organoids wt and 30p from daz 23 diff, 53t from daz 30 diff
% % DB317
% SavePath = 'S:\HCS_Platform\Data\JenniferModamio\2019DB317DB322_analysis\Synuclein\ScriptSyn20190709_JM20190316_DB317_slideNum5_synFil488_syn2568_TH647'; % organoids wt and 30p from daz 23 diff, 53t from daz 30 diff

% DB322 slide 12.1
% SavePath = 'S:\HCS_Platform\Data\JenniferModamio\2019DB317DB322_analysis\Synuclein\ScriptSyn20190709_JM20190402_DB322_slideNum12_SynFilamnt488_SynTotal568_TH647_25_5'; % organoids wt and 30p from daz 23 diff, 53t from daz 30 diff
% DB322 slide 12.2
%SavePath = 'S:\HCS_Platform\Data\JenniferModamio\2019DB317DB322_analysis\Synuclein\ScriptSyn20190709_JM20190402_DB322_slideNum12_SynFilamnt488_SynTotal568_TH647_40_6_testHighSynuclein4level'; % organoids wt and 30p from daz 23 diff, 53t from daz 30 diff

%% batch test old organoids 

% % DB317
% SavePath = 'S:\HCS_Platform\Data\JenniferModamio\2019DB317DB322_analysis\Synuclein\ScriptSyn20190709_JM20190410_DB317_synuclein_322_synuclein_p3_70d_testHighSynuclein4level'; % organoids wt and 30p from daz 23 diff, 53t from daz 30 diff

%% batch from 27/05/19 approx 

% DB322 & 317 30d old 
% SavePath = 'S:\HCS_Platform\Data\JenniferModamio\2019DB317DB322_analysis\Synuclein\ScriptSyn20190709_JM20190529_DB322_DB317_p10p12p14_30dDiff_synFilamnt488_synTotal568_TH647'; % organoids wt and 30p from daz 23 diff, 53t from daz 30 diff

% DB322 & 317 50d old 
%SavePath = 'S:\HCS_Platform\Data\JenniferModamio\2019DB317DB322_analysis\Synuclein\ScriptSyn20190709_JM20190528_DB322_DB317_p8p10p11p12_50dDiff_synFilamnt488_synTotal568_TH647'; % organoids wt and 30p from daz 23 diff, 53t from daz 30 diff

% DB322 & 317 70d old 
% SavePath = 'S:\HCS_Platform\Data\JenniferModamio\2019DB317DB322_analysis\Synuclein\ScriptSyn20190709_JM20190528_DB322_DB317_p6p10p11p12_70dDiff_synFilamnt488_synTotal568_TH647'; % organoids wt and 30p from daz 23 diff, 53t from daz 30 diff

%% batch from 27/05/19 approx 
SavePath = 'S:\HCS_Platform\Data\JenniferModamio\2019DB317DB322_analysis\Synuclein\ScriptSyn20190709_JM20190727_322_317_p10_p11_p12_90dDiff_filmntsyn_488_TotalSyn_56'; % organoids wt and 30p from daz 23 diff, 53t from daz 30 diff



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
%% first batch of organoids 

%%%%%    REMEMEBER CHANGE OBJECTS NAME %%%%%
   
% DB322
% InfoTable = readtable('S:\Operetta\OperettaDB_LCSB\JM20190316_DB322_slideNum6_synFil488_syn2568_TH647\7f484e24-55ec-413e-82d8-e87803c35e54\metadata.csv');
% DB317
% InfoTable = readtable('S:\Operetta\OperettaDB_LCSB\JM20190316_DB317_slideNum5_synFil488_syn2568_TH647_repet2\b1b98016-44d9-4163-901b-3d089f8a73b8\metadata.csv');

% DB322 slide 12.1
% InfoTable = readtable('S:\Operetta\OperettaDB_LCSB\JM20190402_DB322_slideNum12_SynFilamnt488_SynTotal568_TH647_25_5\f9029d9f-dfaa-4139-a23b-18c8e0b841b8\metadata.csv');
% DB322 slide 12.2
%InfoTable = readtable('S:\Operetta\OperettaDB_LCSB\JM20190402_DB322_slideNum12_SynFilamnt488_SynTotal568_TH647_40_6\3caa4a5d-1b2d-4304-9c13-94a9dc0ef2c3\metadata.csv');
% DB317 slide 11
% InfoTable = readtable('S:\Operetta\OperettaDB_LCSB\JM20190402_DB317_slideNum11_synFil488_syn2568_TH647_26_53um_all\87d375ee-dde3-4772-be02-77ad1dc0c5d2\metadata.csv');

%% first batch of organoids 
 % DB322 & 317 70d 
 %  InfoTable = readtable('S:\Operetta\OperettaDB_LCSB\JM20190410_DB317_th_syn_stem_322_th_syn_stem_p3_70d\c2d6ee97-bd40-4a11-b9cf-928440f02090\metadata.csv');

%% batch from 27/05/2019
 
 % DB322 & 317 30d 
 % InfoTable = readtable('S:\Operetta\OperettaDB_LCSB\JM20190604_DB322_DB317_p10p12p14_30dDiff_TH647_TotalSyn568_FilSy\fe24bee8-5b78-465c-b285-754eee662049\metadata.csv');
 % DB322 & 317 50d 
 % InfoTable = readtable('S:\Operetta\OperettaDB_LCSB\JM20190603_DB322_DB317_p8p10p11p12_50dDiff_TH647_TotalSyn568_Fil\b1371aab-d102-497c-b0c6-eb72c0405192\metadata.csv');
 % DB322 & 317 70d 
 % InfoTable = readtable('S:\Operetta\OperettaDB_LCSB\JM20190603_DB322_DB317_p6p10p11p12_70dDiff_TH647_TotalSyn568_Fil\b2093951-4b5e-4a02-b472-85d4357c7518\metadata.csv');

%% batch from 27/07/2019
InfoTable = readtable('S:\Operetta\OperettaDB_LCSB\JM20190727_322_317_p10_p11_p12_90dDiff_filmntsyn_488_TotalSyn_56\31b0a985-085b-4344-be5e-55f5b190368d\metadata.csv');

 
 
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
            %  ObjectsThisOrganoid = f_ImageAnalysisPerOperettaOrganoid_cell_count(Label, CroppedMosaic{1}, CroppedMosaic{2}, CroppedMosaic{3}, CroppedMosaic{4}, ChannelNames, PreviewSavePath);
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
    writetable(Objects, [SavePath, filesep, 'JM20190402_DB322_slideNum12_SynFilamnt488_SynTotal568_TH647_40_6.csv']) %change to name slide 
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




