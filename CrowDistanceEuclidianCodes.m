function  FullSummedMatrix = CodesDistanceEuclidian(LCMIN,PopIN)
warning off
if matlabpool('size') ~= 0
    matlabpool close
end
if exist('CountingDir') ~= 0
    cd('CountingDir')
    delete('f*')
    cd ..
    rmdir('CountingDir')
end
mkdir('CountingDir')

tStart = tic;
Size(1) = size(LCMIN,1);
Size(2) = size(LCMIN,2);
FarmPresIn = zeros(Size(1),Size(2));
FarmPresIn(:,:) = NaN;
% Determine "Farmland" target cells
for Row = 1:Size(1)
    for  Col = 1:Size(2)
        if LCMIN(Row,Col) == 3 ||LCMIN(Row,Col) == 4
            FarmPresIn(Row,Col) = 1;
        end
    end
end
%Prepare distance calculations
display('Caluclating Double matrix')
DoubleSize = Size.*2;
VectorMatrix = zeros((DoubleSize(1).*DoubleSize(2)),2);
for Row = 1:(DoubleSize(1))
    xstart = ((Row-1).*DoubleSize(2))+1;
    xend = Row.*DoubleSize(2);
    VectorMatrix(xstart:xend,1) = Row;
    VectorMatrix(xstart:xend,2) = [1:DoubleSize(2)];
end
DistanceMatrixLarge = PerPointfun(Size(1),Size(2),VectorMatrix,DoubleSize);
save('DoubleMatrix250m','DistanceMatrixLarge')
FullSummedMatrix = zeros(size(LCMIN,1),size(LCMIN,2));
FullSummedMatrix(LCMIN<0) = LCMIN(LCMIN<0);
PopIN(PopIN<0) =NaN;
save('allTestBackup')
clear LCMIN Row Col xstart xend VectorMatrix DoubleSize

%% as parfor for Row
todoRows = Size(1);
clc
[RowArray,EmptyRows] = EmptyFunc(todoRows,FarmPresIn);
matlabpool open Full
parfor Row = 1:todoRows
    List = cell2mat(RowArray(Row));
    if isempty(List) ~= 1
        tic
        LineMatrix = zeros(1,Size(2));            
        Lister = zeros(1,length(List));
        for i = 1:length(List)
            Col = List(i);
            Lister(i)= ExtractPerPointfun(DistanceMatrixLarge,Row,Col,PopIN,Size);
        end
        LineMatrix(1,List) = Lister;
        FullSummedMatrix(Row,:) = FullSummedMatrix(Row,:) + LineMatrix;
        fid = fopen(['CountingDir/f' num2str(Row)],'w');
        fclose(fid);
        Counter = length(dir('CountingDir'))-2;
        fprintf('Finished %i rows with %i cells of a total of %i non empty rows with %3.4f cells per second\n',...
            Counter,length(List),(todoRows-EmptyRows),(((length(List)./toc)).*10))
    end
end
matlabpool close
format short g
tRuns = toc(tStart) 

%%
display('Saving data: Mat-File')
save('FullSummedMatrixEuclidian','FullSummedMatrix')
display('Saving data: Asci-File')
%% 250m
AsciHeading(1,1) = {'ncols         2854'};
AsciHeading(2,1) = {'nrows         5200'};
AsciHeading(3,1) = {'xllcorner     -13520.119510958'};
AsciHeading(4,1) = {'yllcorner     -8.4835694874637'};
AsciHeading(5,1) = {'cellsize      250'};
AsciHeading(6,1) = {'NODATA_value  -9999'};
%% 1KM
% AsciHeading(1,1) = {'ncols         714'};
% AsciHeading(2,1) = {'nrows         1300'};
% AsciHeading(3,1) = {'xllcorner     -13520.119510958'};
% AsciHeading(4,1) = {'yllcorner     -8.4835694874637'};
% AsciHeading(5,1) = {'cellsize      1000'};
% AsciHeading(6,1) = {'NODATA_value  -9999'};

filename = 'myFileEuclidian.asc';
fid = fopen(filename, 'wt');
fprintf(fid, '%s\n%s\n%s\n%s\n%s\n%s\n', char(AsciHeading(1)),char(AsciHeading(2)),...
    char(AsciHeading(3)),char(AsciHeading(4)),char(AsciHeading(5)),...
    char(AsciHeading(6)));  % header
fclose(fid);
dlmwrite(filename,FullSummedMatrix,'delimiter','\t','precision',['%15.',num2str(0),'f'],'-append');
format short g
tEnd = toc(tStart)
end

%%
%One time distance function to twice as large size
function DistanceWeight = PerPointfun(RowIn,ColIn,VectorMatrix,Size)
distances = sqrt((RowIn-VectorMatrix(:,1)).^2 + (ColIn- VectorMatrix(:,2)).^2);
Xin = reshape(distances,Size(2),Size(1))';
modelFun = @(Xin) exp(-(Xin./10)); 
% the function to be used, can be changed accordingly
DistanceWeight = modelFun(Xin);
end
%%
% Shift distance window to applicable
function OutPoint = ExtractPerPointfun(DistanceMatrixLarge,RowIn,ColIn,PopIn,SizeIn)
Start = SizeIn(1)-RowIn+1;
StartCol = SizeIn(2)-ColIn+1;
DistanceWeight = DistanceMatrixLarge(Start:(Start+(SizeIn(1)-1)),...
                    StartCol:(StartCol+(SizeIn(2)-1)));
DistanceWeight(isnan(PopIn)==1) = NaN;
OutPoint = round(nansum(nansum(DistanceWeight.*PopIn)));
end
%%
% Determine which cells in which rows are Farmland; and number of empty
% rows
function [RowArray,EmptyRows] = EmptyFunc(todoRows,FarmPresIn)
FullRows = zeros(todoRows,1);
RowArray = cell(todoRows,1);
for Row = 1:todoRows
    RowArray{Row} = find(FarmPresIn(Row,:)==1);
    if isempty(RowArray{Row}) ~= 1
     FullRows(Row,1) = 1;
    end
end
EmptyRows = todoRows - sum(FullRows);
end

        
%%         
% function DistanceMatrixPoint = DisMatFunc(RowIn,ColIn,~,VectorMatrix,Size)
% DistanceMatrixPoint = zeros(Size(1),Size(2));
% %distances = sqrt((RowIn-VectorMatrixRow).^2 + (bsxfun(@minus,VectorMatrixCol,ColIn)).^2);
% 
% % Return back to matrix format
% 
% % for Listnr = 1:length(distances)
% %     DistanceMatrixPoint(VectorMatrix(Listnr,1),VectorMatrix(Listnr,2)) = distances(Listnr);
% % end
% end


%% as parfor
% count = 0;
% matlabpool open Full
% todoRows = size(LCMIN,1);
% for Row = 1:todoRows
%     List = find(FarmPresIn(Row,:)==1);
%     if isempty(List) ~= 1
%         LineMatrix = zeros(1,size(LCMIN,2));            
%         count = count + length(List);
%         clc
%         fprintf('%i th total count for row %i of %i Farmland cells with %f cells per second\n', count, Row, length(List),(count./toc))
%         Lister = zeros(1,length(List));
%         parfor i = 1:length(List)
%             Col = List(i);
%             Lister(i)= PerPointfun(LCMIN,Row,Col,PopIN,VectorMatrix,Size);
%         end
%         LineMatrix(1,List) = Lister;
%         clear Lister
%         FullSummedMatrix(Row,:) = FullSummedMatrix(Row,:) + LineMatrix;
%     end
% end
% matlabpool close

%% as job
%count = 1;
% for Row = 1:400%size(LCMIN,1)
%     List = find(FarmPresIn(Row,:)==1);
%     if isempty(List) ~= 1
%         job = createJob('configuration','Full');
%         for i = 1:length(List)
%             if i == length(i)
%                 clc
%                 fprintf('%i th total count for row %i of %i Farmland cells \n', count, Row, length(List))
%             end
%             Col = List(i);
%             createTask(job, @PerPointfun, 1,{LCMIN,Row,Col,PopIN});
%             count = count + 1;
%         end
%         submit(job);
%         waitForState(job, 'finished');
%         Results = getAllOutputArguments(job);
%         for i = 1:length(List)
%             Col = List(i);
%             FullSummedMatrix(Row,Col) = FullSummedMatrix(Row,Col) + cell2mat(Results(i,1));
%         end
%         destroy(job)
%     end
%end

%% 
% 
% display('Running Prediction statistics')
% display('Running loops')
% stints_total = ceil(List.loop_max_Predict./100);
% % so each number not dividedable by ten will be scaled
% % upwards, for that the actual done loops are counted (see
% % below and used as max loops value below
% for stint = 1:1: stints_total
%     job = createJob('configuration','Full');
%     for withinloop = 1:100
%         loop = ((stint-1).*100)+ withinloop;
%         clc
%         display('Running Predictive runs')
%         display(List.Outputs(trait))
%         display(traittext)
%         display(loop)
%         createTask(job, @PredictLoop, 0,{loop,List,PropRadi,trait,traittype});
%         %Loops. Means all transfer parameters to keep should be in
%         %the seperate swap files
%     end
%     submit(job);
%     time_out_time = 600;
%     waitForState(job, 'finished',time_out_time);
%     destroy(job)
% end

