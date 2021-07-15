%% JM 06.2020

%% SNCA lines DB322 DB317 

%% User inputs
clear


SetupMode = 0; % 1 for creating numeric organoid labels OR 0 for linking the final analysis to human labels


% 23052020_30d_hMO_number1__condition6_MAP2_647_TH_488_TUJ1_568_p9_10_11_JM_exp29
%SlideLayout = 'exp29.txt'; % CAREFULL! old acquisition. adjust thresholds 
% 24052020_70d_hMO_number2__condition6_MAP2_647_TH_488_TUJ1_568_p9_10_11_JM_exp30
%SlideLayout = 'exp30.txt'; 
%  06062020_90d_hMO_number1__condition6_MAP2_647_TH_488_TUJ1_568_p9_10_11_JM_exp40
SlideLayout = 'exp40.txt'; 
%06062020_90d_hMO_number1_exp123_condition6_MAP2_647_TH_488_TUJ1_568_p9_10_11_JM_exp48
%SlideLayout = 'exp48.txt'; 
% 25062020_70d_hMO_number1_exp123_condition6_MAP2_647_TH_488_TUJ1_568_p9_10_11_JM_exp49
%SlideLayout = 'exp49.txt';
% 25062020_30d_hMO_number1_exp123_condition6_MAP2_647_TH_488_TUJ1_568_p9_10_11_JM_exp50
%SlideLayout = 'exp50.txt'; 
%25062020_50d_hMO_number1_exp123_condition6_MAP2_647_TH_488_TUJ1_568_p9_10_11_JM_exp51
%SlideLayout = 'exp51.txt'; 
% 25062020_15d_hMO_number1_exp123_condition6_MAP2_647_TH_488_TUJ1_568_p9_10_11_JM_exp52
%SlideLayout = 'exp52.txt'; 

% 13072020_70dp1_hMO_number1_exp567_condition6_MAP2_647_TH_488_TUJ1_568_p9_10_11_JM_exp66
%SlideLayout = 'exp66.txt';
%13072020_70dp2_hMO_number1_exp567_condition6_MAP2_647_TH_488_TUJ1_568_p9_10_11_JM_exp67
%SlideLayout = 'exp67.txt';
%14072020_50dp1_hMO_number1_exp567_condition6_MAP2_647_TH_488_TUJ1_568_p9_10_11_JM_exp69
%SlideLayout = 'exp69.txt';

% 18072020_50dp2_hMO_number1_exp567_condition6_MAP2_647_TH_488_TUJ1_568_p9_10_11_JM_exp70
%SlideLayout = 'exp70.txt';
% 18072020_150-180d_hMO_number1_exp123_condition6_MAP2_647_TH_488_TUJ1_568_p9_10_11_JM_exp71
%SlideLayout = 'exp71.txt';
%18072020_15d_hMO_number1_exp567_condition6_MAP2_647_TH_488_TUJ1_568_p9_10_11_JM_exp72
%SlideLayout = 'exp72.txt';
%13072020_30dp1_hMO_number1_exp567_condition6_MAP2_647_TH_488_TUJ1_568_p9_10_11_JM_exp73
%SlideLayout = 'exp73.txt';
%13072020_90dp2_hMO_number1_exp567_condition6_MAP2_647_TH_488_TUJ1_568_p9_10_11_JM_exp74
%SlideLayout = 'exp74.txt';

%SlideLayout = 'exp75.txt';


%% 
%SavePath = 'S:\HCS_Platform\Data\JenniferModamio\2020\hMO_characterisation\JM20200630_SNCA_Condition6_TUJ1_MAP2_TH\JM_exp40_2_'; % organoids wt and 30p from daz 23 diff, 53t from daz 30 diff

%SavePath = 'S:\HCS_Platform\Data\JenniferModamio\2020\hMO_characterisation\JM20200630_SNCA_Condition6_TUJ1_MAP2_TH\JM_exp29_3'; % organoids wt and 30p from daz 23 diff, 53t from daz 30 diff

%SavePath = 'S:\HCS_Platform\Data\JenniferModamio\2020\hMO_characterisation\JM20200630_SNCA_Condition6_TUJ1_MAP2_TH\JM_exp30_2_'; 


%SavePath = 'S:\HCS_Platform\Data\JenniferModamio\2020\hMO_characterisation\JM20200630_SNCA_Condition6_TUJ1_MAP2_TH\JM_exp48_3'; % organoids wt and 30p from daz 23 diff, 53t from daz 30 diff
%SavePath = 'S:\HCS_Platform\Data\JenniferModamio\2020\hMO_characterisation\JM20200630_SNCA_Condition6_TUJ1_MAP2_TH\JM_exp51_2_3_3'; 
SavePath = 'S:\HCS_Platform\Data\JenniferModamio\2020\hMO_characterisation\JM20200630_SNCA_Condition6_TUJ1_MAP2_TH\JM_exp40_2021_test'; 
%SavePath = 'S:\HCS_Platform\Data\JenniferModamio\2020\hMO_characterisation\JM20200630_SNCA_Condition6_TUJ1_MAP2_TH\JM_exp48'; 
%SavePath = 'S:\HCS_Platform\Data\JenniferModamio\2020\hMO_characterisation\JM20200630_SNCA_Condition6_TUJ1_MAP2_TH\JM_exp30'; 
%SavePath = 'S:\HCS_Platform\Data\JenniferModamio\2020\hMO_characterisation\JM20200630_SNCA_Condition6_TUJ1_MAP2_TH\JM_exp52_2_Example'; 
%SavePath = 'S:\HCS_Platform\Data\JenniferModamio\2020\hMO_characterisation\JM20200630_SNCA_Condition6_TUJ1_MAP2_TH\JM_exp50_2_'; 
%SavePath = 'S:\HCS_Platform\Data\JenniferModamio\2020\hMO_characterisation\JM20200630_SNCA_Condition6_TUJ1_MAP2_TH\JM_exp49_3'; 
%SavePath = 'S:\HCS_Platform\Data\JenniferModamio\2020\hMO_characterisation\JM20200630_SNCA_Condition6_TUJ1_MAP2_TH\JM_exp70_3'; 
%SavePath = 'S:\HCS_Platform\Data\JenniferModamio\2020\hMO_characterisation\JM20200630_SNCA_Condition6_TUJ1_MAP2_TH\JM_exp73_3'; 
%SavePath = 'S:\HCS_Platform\Data\JenniferModamio\2020\hMO_characterisation\JM20200630_SNCA_Condition6_TUJ1_MAP2_TH\JM_exp74_3'; 
%SavePath = 'S:\HCS_Platform\Data\JenniferModamio\2020\hMO_characterisation\JM20200630_SNCA_Condition6_TUJ1_MAP2_TH\JM_exp72_3_7on'; 


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
    
%InfoTable = readtable('S:\Operetta\OperettaDB_LCSB\JM_exp40\6a6fda86-9e45-4a9b-9651-a947e6766d0b\metadata.csv');
%InfoTable = readtable('S:\Operetta\OperettaDB_LCSB\JM_exp29\a877a635-b6ef-45e2-b335-b0534ca46134\metadata.csv');
%InfoTable = readtable('S:\Operetta\OperettaDB_LCSB\JM_exp30\7c594bbe-463b-4dfd-873a-175f4ba15fd9\metadata.csv');

%InfoTable = readtable('S:\Operetta\OperettaDB_LCSB\JM_exp29_2\335e267b-dcfb-4bb4-96d0-626ae7b14731\metadata.csv');
InfoTable = readtable('S:\Operetta\OperettaDB_LCSB\JM_exp40_2\62de1695-dd9c-4051-a3a3-5f714f0b4c33\metadata.csv');
%InfoTable = readtable('S:\Operetta\OperettaDB_LCSB\JM_exp30_2\6d98b730-70a8-427f-8242-fcd18cd1c270\metadata.csv');

%InfoTable = readtable('S:\Operetta\OperettaDB_LCSB\JM_exp51\d82af3a2-0644-4fea-a070-c0270fdd8b31\metadata.csv')
%InfoTable = readtable('S:\Operetta\OperettaDB_LCSB\JM_exp40\6a6fda86-9e45-4a9b-9651-a947e6766d0b\metadata.csv');
%InfoTable = readtable('S:\Operetta\OperettaDB_LCSB\JM_exp48\705859e7-5127-45e1-916c-046f2844d48a\metadata.csv');
% InfoTable = readtable('S:\Operetta\OperettaDB_LCSB\JM_exp30\7c594bbe-463b-4dfd-873a-175f4ba15fd9\metadata.csv');
% InfoTable = readtable('S:\Operetta\OperettaDB_LCSB\JM_exp52\7975e381-a852-47b4-9383-20513d0d1ae3\metadata.csv');
%InfoTable = readtable('S:\Operetta\OperettaDB_LCSB\JM_exp50\71fe4cf8-c4bb-47a7-a9ba-44b7a2179d3a\metadata.csv');
%InfoTable = readtable('S:\Operetta\OperettaDB_LCSB\JM_exp49\42c96d71-30af-4958-81b4-3d2aab83da57\metadata.csv');

%InfoTable = readtable('S:\Operetta\OperettaDB_LCSB\JM_exp70\0e7d6971-173a-4642-82df-1720be9a6b64\metadata.csv');
%InfoTable = readtable('S:\Operetta\OperettaDB_LCSB\JM_exp73\3ebf42b9-0c51-4ea8-8214-6f9becac5f95\metadata.csv');
%InfoTable = readtable('S:\Operetta\OperettaDB_LCSB\JM_exp72\afbc5747-be89-4acd-87ae-50a91f54a494\metadata.csv');



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
   
   Groups =6
    for g = Groups % Number of organoids
    g=6
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
                        ImPathThisPlane = strrep(ImPathThisPlane, '\\atlas.uni.lux\LCSB_HCS\Operetta', 'S:\Operetta');
                        
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
            ObjectsThisOrganoid = f_ImageAnalysisPerOperettaOrganoid_neuronalRatio(Label, CroppedMosaic{1}, CroppedMosaic{2}, CroppedMosaic{3}, CroppedMosaic{4}, ChannelNames, PreviewSavePath);
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

    for g =  Groups
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




