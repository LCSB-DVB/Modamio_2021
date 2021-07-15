%% JM 13.03.2019

%% SNCA lines 

%% stainings: 
% 421 total synuclein ms 1:1000 568
% 439 syn filament rb 1:5000 488
% 364 MAP2 1:1000 647 

%% User inputs
clear
clc

SetupMode = 0; % 1 for creating numeric organoid labels OR 0 for linking the final analysis to human labels

% SlideLayout = '17022021_50d_hMO_exp12_c1_synFil_total_right.txt'; 
% SavePath = 'S:\HCS_Platform\Data\JenniferModamio\2021\JM20210223_SNCA_GDNF_C1_SYN_SynF_MAP2\17022021_50d_hMO_exp12_c1_synFil_total_right'; % organoids wt and 30p from daz 23 diff, 53t from daz 30 diff

 SlideLayout = '17022021_50d_hMO_exp12_c1_synFil_total_left.txt'; 
 SavePath = 'S:\HCS_Platform\Data\JenniferModamio\2021\JM20210223_SNCA_GDNF_C1_SYN_SynF_MAP2\JM_17022021_50d_hMO_exp12_c1_synFil_total_left_mask'; % organoids wt and 30p from daz 23 diff, 53t from daz 30 diff

% SlideLayout = 'JM_17022021_50d_hMO_exp11_c1_synPhosp_total.txt'; 
% SavePath = 'S:\HCS_Platform\Data\JenniferModamio\2021\JM20210223_SNCA_GDNF_C1_SYN_SynF_MAP2\JM_17022021_50d_hMO_exp11_c1_synPhosp_total'; % organoids wt and 30p from daz 23 diff, 53t from daz 30 diff

%  SlideLayout = 'JM_17022021_50d_hMO_exp10_c1_synPhosp_total_GDNFp_neg.txt'; 
%  SavePath = 'S:\HCS_Platform\Data\JenniferModamio\2021\JM20210223_SNCA_GDNF_C1_SYN_SynF_MAP2\JM_17022021_50d_hMO_exp10_c1_synPhosp_total_GDNFp_neg'; % organoids wt and 30p from daz 23 diff, 53t from daz 30 diff
% 
%  SlideLayout = 'JM_17022021_50d_hMO_exp10_c1_control.txt'; 
%  SavePath = 'S:\HCS_Platform\Data\JenniferModamio\2021\JM20210223_SNCA_GDNF_C1_SYN_SynF_MAP2\JM_17022021_50d_hMO_exp10_c1_synPhosp_total_control_masks'; % organoids wt and 30p from daz 23 diff, 53t from daz 30 diff


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

% InfoTable = readtable('S:\Operetta\OperettaDB_LCSB\17022021_50d_hMO_exp12__condition1_synFil_488_SynTotal_568_right\4f2ce1b1-19a8-49cd-8d41-9741589ca55c\metadata.csv');
 InfoTable = readtable('S:\Operetta\OperettaDB_LCSB\17022021_50d_hMO_exp12__condition1_synFil_488_SynTotal_568_left\c8ba3ae5-fe5c-431e-a5fd-7e64515d58a3\metadata.csv');

%InfoTable = readtable('S:\Operetta\OperettaDB_LCSB\JM_17022021_50d_hMO_exp11__condition1_synFil_488_SynTotal_568_MA\42d31fd0-2997-4180-8d15-72e06cebed99\metadata.csv');
% InfoTable = readtable('S:\Operetta\OperettaDB_LCSB\JM_20032021_50d_hMO_exp10_c1_synFil_488_SynTotal_568_MAP2_647_co\9150314d-05ff-4a53-91ce-f5a015199b2d\metadata.csv');

% InfoTable = readtable('S:\Operetta\OperettaDB_LCSB\12032021_50d_hMO_exp10_condition1_synFil_488_SynTotal_568_MAP2_6\c7483031-0c22-4d24-922c-1a3b25bec215\metadata.csv');

 
 
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
            ObjectsThisOrganoid = f_ImageAnalysisPerOperettaOrganoid_Synuclein(Label, CroppedMosaic{1}, CroppedMosaic{2}, CroppedMosaic{3}, CroppedMosaic{4}, ChannelNames, PreviewSavePath);
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
    
    % there was onlz 5 objects. an error was encounter in the second
    % organoid (corresponding to A53T?own derivation). I re'run the loop
    % with g being 2. I save it in this separeted folder to do not
    % overwrite the previous results and have a back'up in case is needed. 
    
    
%     Objects = vertcat(ObjectsAll{:});
%     save([SavePath, filesep, 'Objects2.mat'], 'Objects');
%     writetable(Objects, [SavePath, filesep, 'Objects2.csv'])
%     writetable(Objects, [SavePath, filesep, 'Objects2.xlsx'])
%     % saved objectsAll, saved error

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




