% This is the population distance decay gravity  model based on LCM 
% categories,% population density and potentoial traveling time trough 
% cells. Used is the Schlaepfer distance decay function (line 180).
% Inputs are a 1) Population density raster; 2) an LCM raster; 3) a speed
% raster based on road presence (at a larger resolution; 4) a narrow
% resolution distance to main road raster.

%Main module with NumIn is resolution that is matched to input data
% FrequencyIn is the frequency of visit (1=Yearly; 12=Monthly; 52=Weekly)
function  FullSummedMatrix = CodesRoadDistance(NumIn,StartTranch,FrequencyIn)
warning off

% ********************
% THE START
fprintf('Starting preparations \n')
if matlabpool('size') ~= 0
    matlabpool close
end
Transfer.NumIn = NumIn;
Transfer.StartTranch = StartTranch;
Transfer.Frequency = FrequencyIn;

if Transfer.Frequency == 1
    Transfer.FrequencyText = 'Yearly';
elseif Transfer.Frequency == 12
    Transfer.FrequencyText = 'Monthly';
elseif Transfer.Frequency == 52
    Transfer.FrequencyText = 'Weekly';
else
    Transfer.FrequencyText = mat2str(datenum(Transfer.Frequency));
end
clear NumIn StartTranch FrequencyIn

if exist('FullSummedMatrixTmp.mat') ~= 0 && Transfer.StartTranch == 1
    [Transfer.StartTranch] = TemporaryFileAnswer;
    if Transfer.StartTranch == 1
        delete('FullSummedMatrixTmp.mat'); % delete the temporary file
    end
end
% Initiate counting dir (so progress can be counted in a parfor loop)
if exist('CountingDir') ~= 0 && Transfer.StartTranch == 1
    cd('CountingDir')
    delete('f*')
    cd ..
    rmdir('CountingDir')
end
mkdir('CountingDir')

%% Load applicable data and set asci (-9999) to NaN
[LCMIN,PopIN,AttractiveRaster,Transfer,DistanceRastExtra,AsciHeading] =...
    LoadData(Transfer);

%% Define Farmland
fprintf('\n Defining Farmland \n')
Transfer.Size(1) = size(LCMIN,1);
Transfer.Size(2) = size(LCMIN,2);
FarmPresIn = zeros(Transfer.Size(1),Transfer.Size(2));
FarmPresIn(:,:) = NaN;
FarmPresIn(LCMIN >0) = 1; 

% Kick out certain categories if needed
% FarmPresIn(LCMIN<20) = NaN;

% Set non existing land use to NaN;
FarmPresIn(isnan(LCMIN)==1) = NaN;
%% Print total number of cells
Transfer.NumberOfCells = nansum(nansum(FarmPresIn));
fprintf('\n Total number of cells to run: %i \n', Transfer.NumberOfCells);
%% Redefine attractiveness raster
AttractiveRaster(AttractiveRaster== 10) = NaN; % No entry strict nature reserves
AttractiveRaster = AttractiveRaster./5; % make it a proportion
%%
fprintf('\n Define full and empty cells \n')
Transfer.todoRows = Transfer.Size(1);
[Transfer.RowArray,Transfer.EmptyRows,Transfer.NrCellsRow] = EmptyFunc(Transfer.todoRows,FarmPresIn); 
save('LCDATA','LCMIN', 'FarmPresIn') % save for later use and clear
clear FarmPresIn TwoFiveKList DistanceRast250mExtra LCMIN DistanceRast1kmExtra
%% Main run as parfor for Row within tranches
% MAIN RUN START
%load potential earlier file
if Transfer.StartTranch == 1
    FullSummedMatrix = zeros(Transfer.Size(1),Transfer.Size(2));
    Transfer.TimeSpendEarlierRun = 0;
else
    load('FullSummedMatrixTmp.mat')
    Transfer.TimeSpendEarlierRun = TimeSpend; %#ok<NODEF>
end
matlabpool open 8
Transfer.ProcessorstoUse = matlabpool('size'); 

% Run actual calculations as double loop. The double loop is for security
% to save inbetween data. Number of tranches is pure arbitrary
% Determine nr rows per tranch 
fprintf('\n Starting per cell run \n')
Transfer.NumRow = ceil(Transfer.todoRows./Transfer.NumTranch);
Transfer.PercDone = 0;
Transfer.End = 0;
Transfer.tStart = tic;
save('allTestBackup') % for testing input data
for trh = Transfer.StartTranch:Transfer.NumTranch
    % determine the rows per tranch
    Transfer = SetRowsToRun(Transfer,trh);
    % Run as parfor loop within tranches, the output file is filled up accordingly
    FullSummedMatrix = RunAllCells(FullSummedMatrix,Transfer,PopIN,DistanceRastExtra,...
        AttractiveRaster);
    TimeSpend = toc(Transfer.tStart); %#ok<NASGU>
    % save the inbetween data
    save('FullSummedMatrixTmp','FullSummedMatrix','TimeSpend');
end
if matlabpool('size') ~= 0
    matlabpool close
end
% MAIN RUN END

%% Make sure only true cells are shown
load('LCDATA.mat')
delete('LCDATA.mat');
%FullSummedMatrix(LCMIN<0) = LCMIN(LCMIN<0);
FullSummedMatrix(isnan(FarmPresIn)==1) = NaN;
clear FarmPresIn LCMIN

% Print time used
Transfer.End = 1;
[~,timevar] = TimeLeftShow(Transfer);
fprintf('\n This has run has taken %i hours, %i minutes and %i seconds \n',timevar.hr,timevar.mn,timevar.sec)

%% Saving all data
fprintf('\n Saving data: Mat-File \n')
name = ['FullSummedMatrix_',mat2str(Transfer.NumIn),'_',Transfer.FrequencyText];
save(name,'FullSummedMatrix')
delete('FullSummedMatrixTmp.mat'); % delete the temporary file
delete('allTestBackup.mat'); % delete the temporary file

% asci file save with heading
filename = ['PressureRaster',mat2str(Transfer.NumIn),'_',Transfer.FrequencyText,mat2str(datenum(date)),'.asc'];
SaveAsci(filename,FullSummedMatrix,AsciHeading,Transfer.AsciRound)

% Print total time used
[~,timevar] = TimeLeftShow(Transfer);
fprintf('\n This full time was %i hours, %i minutes and %i seconds \n',timevar.hr,timevar.mn,timevar.sec)
end
% THE END
% ***************************


%% SUBFUNCTIONs
%% Main Function
function FullSummedMatrix = RunAllCells(FullSummedMatrix,Transfer,PopIN,DistanceRastExtra,...
                 AttractiveRaster)                        
fprintf('\n Restarting with tranch number %i with row number %i of total %i \n',Transfer.tranch,Transfer.RowStart, Transfer.todoRows)
parfor Row = Transfer.RowStart:Transfer.RowEnd
    List = cell2mat(Transfer.RowArray(Row)); %#ok<PFBNS>
    if isempty(List) ~= 1
        tic
        LineMatrix = zeros(1,Transfer.Size(2));             
        Lister = zeros(1,length(List));
        DHIDStore = -1;
        StoreWeight = [];
        for i = 1:length(List)
            Col = List(i);
            [Lister(i),StoreWeight,DHIDStore] =...
                ExtractPerPointfunRoads(Row,Col,PopIN,Transfer.RowColList,DistanceRastExtra,...
                 AttractiveRaster,DHIDStore,StoreWeight,Transfer.NumIn,Transfer.Frequency);
        end
        LineMatrix(1,List) = Lister;
        FullSummedMatrix(Row,:) = FullSummedMatrix(Row,:) + LineMatrix; 
        
        %CountProcedure  
        fid = fopen(['CountingDir/f' num2str(Row)],'w');
        fclose(fid);
        Counter = length(dir('CountingDir'))-2;
        if (Row./5) == ceil(Row./5) %#ok<BDSCA>
            fprintf('Finished %i rows with %i cells of a total of %i non empty rows with %3.4f cells per second\n',...
                Counter,length(List),(Transfer.todoRows-Transfer.EmptyRows),...
                (((length(List)./toc)).*Transfer.ProcessorstoUse))
        end
    end
end
end
%% Calculate gravity distance outcome over all source rasters per target cell
function [OutPoint,StoreWeight,DHIDStore] = ExtractPerPointfunRoads(RowIn,ColIn,...
            PopIn,GridList,DistanceRastIn,AttractiveRaster,DHIDStore,StoreWeight,NumIn,Frequency)
if NumIn == 250
    SizeFactor = 10; % for 2_5km
elseif NumIn == 1000
    SizeFactor = 5; % for 5km
end
ToTake = find(GridList.Row == ceil(RowIn./SizeFactor));
if isempty(ToTake) ~= 1
    Test = GridList(ToTake,:);
    ToTake = find(Test.Col == ceil(ColIn./SizeFactor));
    if isempty(ToTake) ~= 1
        DHID =  Test.DH_ID(ToTake);
    end
end
if isempty(ToTake) ~= 1
    if DHID ~= DHIDStore
        if NumIn == 250
            file = ['DistanceRoad2_5kmAsci/distanceroad2_5_',mat2str(DHID),'.asc'];
        elseif NumIn == 1000
            file = ['DistanceRoad5kmAsci/distancepoint_',mat2str(DHID),'.asc'];
        end
        %display(file)
        [DistanceRastLarge] = cut_top_off_asciiInHere(file);
        Sizes(1) =  (ceil((size(DistanceRastIn,1)./SizeFactor))).*SizeFactor;
        Sizes(2) =  (ceil((size(DistanceRastIn,2)./SizeFactor))).*SizeFactor;
        DistanceRaster = zeros(Sizes(1),Sizes(2));
        DistanceRaster(:,:) = NaN;
        for Row = 1:1:size(DistanceRastLarge,1)
            A = DistanceRastLarge(Row,:);
            RowStart = ((Row-1).*SizeFactor)+1;
            RowEnd = Row.*SizeFactor;
            Tmp1 = (A(ceil((1:(numel(A))*SizeFactor)/SizeFactor)));
            DistanceRaster(RowStart:RowEnd,:) = repmat(Tmp1(:).',SizeFactor,1); %#ok<*BDSCI>
        end
        clear DistanceRastLarge Row
        if NumIn == 250
            if DHID <= 10
                DistanceRaster = DistanceRaster./1000;
            end
            % since it was set at meters not km as the rest.
            DistanceRaster(:,2855:2860) = []; %Only for 2.5km,
        %Size needed = (5200,2854); checked it is the last cols that are extra
        %above the LCM grid!, which is correct since Fishnet starts South-west
        %corner based on the LCM and no rows need removing
        elseif NumIn == 1000
            DistanceRaster(:,715) = []; %Only for 5km
        end
        DistanceRaster(DistanceRaster<0) = NaN; % remove -9999 = NaN;
        DistanceRaster = DistanceRaster + DistanceRastIn;
        clear DistanceRastIn
        
        % Schlaepfer distance decay function
        factor = 2.17;
        %modelFun = @(DistanceRaster) AttractiveRaster(RowIn,ColIn)./((Frequency.*(DistanceRaster+1)).^factor);
        % the function to be used, can be changed accordingly
        DistanceWeight = 1./((Frequency.*(DistanceRaster+1)).^factor);
        %DistanceWeight(DistanceWeight>1) = 1;
        DistanceWeight(isnan(PopIn)==1) = NaN;
        %Store parameters
        StoreWeight = DistanceWeight;
        DHIDStore = DHID;
    else
        DistanceWeight = StoreWeight;
    end
    DistanceWeight(RowIn,ColIn) = 1;
    Roundfactor = (1000./NumIn)^2; % Introduce a scaling factor to make it competabile to 1-km;
    % Combine with PopSize and atrractiveness
    OutPoint = (AttractiveRaster(RowIn,ColIn).*(nansum(nansum(DistanceWeight.*PopIn))));
    OutPoint = (round(OutPoint.*Roundfactor))./Roundfactor;
    if isinf(OutPoint) ==1
        display(RowIn)
        display(ColIn)
        save('all')
        ccc
    end
else
    OutPoint = 0;
end
end

%% load data function
function [LCMIN,PopIN,AttractiveRaster,Transfer,DistanceRastExtra,AsciHeading] =...
    LoadData(Transfer)
if Transfer.NumIn == 250 % 250 meter resolution
    load('UK250mRecreationData.mat')
    LCMIN = UKLCM250m;
    PopIN = UKPop250m;
    AttractiveRaster = UKattractiveness;
    %Prepare for distance calculations % 2_5 km
    Transfer.RowColList = TwoFiveKList;
    Transfer.NumTranch = 50;
    Transfer.AsciRound = 2;
    DistanceRastExtra = DistanceRast250mExtra;
    DistanceRastExtra(:,2855:2860) = []; %Only for 2.5km,
    %Size needed = (5200,2854); checked it is the last cols that are extra
    %above the LCM grid!, which is correct since Fishnet starts South-west
    %corner based on the LCM and no rows need removing 
    clear UKLCM250m UKPop250m UKattractiveness DistanceRast250mExtra TwoFiveKList SlopesAbove250m
    %% 250m heading
    AsciHeading(1,1) = {'ncols         2854'};
    AsciHeading(2,1) = {'nrows         5200'};
    AsciHeading(3,1) = {'xllcorner     -13520.119510958'};
    AsciHeading(4,1) = {'yllcorner     -8.4835694874637'};
    AsciHeading(5,1) = {'cellsize      250'};
    AsciHeading(6,1) = {'NODATA_value  -9999'};
elseif Transfer.NumIn == 1000 || Transfer.NumIn == -999 % 1-km resolution is Transfer.NumIn -999 is a test mode
    load('UK1kmRecreationData.mat')
    LCMIN = UKLCM1km;
    PopIN = UKPop1km;
    AttractiveRaster = UKattractiveness;
    %Prepare for distance calculations % 5 km
    Transfer.RowColList = FiveKList;
    Transfer.NumTranch = 20;
    Transfer.AsciRound = 0;
    DistanceRastExtra = DistanceRast1kmExtra;    
    clear UKLCM1km UKPop1km UKattractiveness FiveKList DistanceRast1kmExtra SlopesAbove1km
    %% 1KM heading
    AsciHeading(1,1) = {'ncols         714'};
    AsciHeading(2,1) = {'nrows         1300'};
    AsciHeading(3,1) = {'xllcorner     -13520.119510958'};
    AsciHeading(4,1) = {'yllcorner     -8.4835694874637'};
    AsciHeading(5,1) = {'cellsize      1000'};
    AsciHeading(6,1) = {'NODATA_value  -9999'};
end

%% Set NaN values (instead of -9999)
LCMIN(LCMIN<=0) = NaN; %zero values cannot exist
AttractiveRaster(AttractiveRaster<=0) = NaN; %zero values cannot exist
PopIN(PopIN<0) =NaN; %zero values are true data
DistanceRastExtra(DistanceRastExtra<0) = NaN; %zero values are true data
end

%% Determine which cells in which rows are Farmland; and number of empty rows
function [RowArray,EmptyRows,NrCellsperRow] = EmptyFunc(todoRows,FarmPresIn)
FullRows = zeros(todoRows,1);
NrCellsperRow = zeros(todoRows,1);
RowArray = cell(todoRows,1);
for Row = 1:todoRows
    List = find(FarmPresIn(Row,:)==1);
    RowArray{Row} = List;   

    if isempty(RowArray{Row}) ~= 1
        FullRows(Row,1) = 1;
        NrCellsperRow(Row,1) = length(List);
    end
    clear List
end
EmptyRows = todoRows - sum(FullRows);
end
%% determine the rows per tranch
function Transfer = SetRowsToRun(Transfer,trh)
Transfer.tranch = trh;
Transfer.RowStart = 1+ ((trh-1).*Transfer.NumRow);
Transfer.RowEnd = trh.*Transfer.NumRow;
if Transfer.RowEnd >= Transfer.todoRows
    Transfer.RowEnd = Transfer.todoRows;
end
[Transfer,timevar] = TimeLeftShow(Transfer);
fprintf('\n Done is %3.2f percent of a total of %i cells with %3.2f cells per second',...
    Transfer.PercDone, Transfer.NumberOfCells, (1./timevar.TimePerCell))
fprintf('\n Estimated time left is %i hours and %i minutes  \n',timevar.hr, timevar.mn)
end

%% Cut the heading from an asci file
function [VargOut] = cut_top_off_asciiInHere(NameIn)
warning off
number_of_heading_lines = 6;
fid = fopen(NameIn,'r');
InputText=textscan(fid,'%s',number_of_heading_lines,'delimiter','\n');
for head = 1:1:number_of_heading_lines
    temp = InputText{1,1}{head,1};
    length_head = length(temp);
    temp_text = temp(1,15:length_head);
    Intro(head) = str2num(temp_text);
    %IntroHead{head,1} = temp;
    clear temp
    clear length_head
    clear temp_text
end
%save('IntroHead','IntroHead');
ncols = Intro(1); %x
nrows = Intro(2); %y
Block = 1;   % Initialize block index
%sprintf('Block: %s', num2str(Block));                % Display block number
InputText=textscan(fid,'%f32','delimiter',','); % Read data block
% Convert to numerical array from cell
Data{Block,1}=cell2mat(InputText); %#ok<*SAGROW> 
data_array = Data{1,1};
clear Data

for row = 1:nrows
    start = ((row-1) *ncols)+1;
    ends = ((row-1) *ncols)+ ncols;
    VargOut(row,:) = data_array(start:ends);
end

fclose('all');
clearvars -except VargOut
end

%% Save asci files with pre defined heading
function SaveAsci(filenameIn,MatrixIn,AsciHeading,AsciRound)
% make sure the target asci file doesn't exist anymore
% else it may write badly.
delete(filenameIn) 

MatrixIn(isnan(MatrixIn)==1) = -9999; % for the ascii file different No data numbers
fprintf('\n Saving data: Asci-File \n')
fid = fopen(filenameIn, 'wt');
fprintf(fid, '%s\n%s\n%s\n%s\n%s\n%s\n', char(AsciHeading(1)),...
    char(AsciHeading(2)), char(AsciHeading(3)),char(AsciHeading(4)),...
    char(AsciHeading(5)), char(AsciHeading(6)));  % header
fclose(fid);
dlmwrite(filenameIn,MatrixIn,'delimiter','\t','precision',...
            ['%15.',num2str(AsciRound),'f'],'-append');
end

%% Time show function
function [Transfer,timevar] = TimeLeftShow(Transfer)
if Transfer.tranch > 0
    endPrevious = (Transfer.tranch-1).*Transfer.NumRow;
    CellsDone = (nansum(nansum(Transfer.NrCellsRow(1:endPrevious,1))));
    Transfer.PercDone = ((CellsDone./Transfer.NumberOfCells).*100);
    tRan = toc(Transfer.tStart)+ Transfer.TimeSpendEarlierRun;
    timevar.TimePerCell = tRan./CellsDone;
    if Transfer.End == 0
        EstimatedTimeLeft = (Transfer.NumberOfCells.*timevar.TimePerCell)-tRan;
    else
        EstimatedTimeLeft = tRan;
    end
    timevar.hr=fix((EstimatedTimeLeft/3600));
    timevar.mn = fix(((EstimatedTimeLeft/3600)-timevar.hr).*60);
    timevar.sec=round(EstimatedTimeLeft - ((timevar.hr.*3600) + (timevar.mn.*60)));
end
end
%% Check for existing temporory file so one can start at any point
function [StartTranch] = TemporaryFileAnswer
reply = 'N';
while reply ~= 'Y'
    reply = input('      Are you sure you want to remove the previous temporrary file (Y or N): ', 's');
    if strcmp('y',reply) == 1
        reply = 'Y';
    end
    if strcmp('Y',reply) == 1
        StartTranch = 1;
    elseif strcmp('N',reply) == 1 || strcmp('n',reply) ==1
        replyNr = input('      Enter the tranch number of the starting row or press control-c to abort (number): ');
        while isnumeric(replyNr) ~= 1 || (round(replyNr) ~= replyNr)
            str = sprintf('      False entry, enter integer');
            disp(str)
            replyNr = input('      Enter the tranch number of the starting row or press control-c to abort (number): ');
        end
        reply = 'Y';
        StartTranch = replyNr;
        str = sprintf('\n      This run will start at tranch %i \n',StartTranch);
        disp(str)
    elseif strcmp('N',reply) ~= 1 && strcmp('n',reply) ~=1 && strcmp('Y',reply) ~= 1 && strcmp('y',reply) ~=1
        str = sprintf('      False entry: %s ', reply);
        disp(str)
        reply = 'N';
    end
end
end

