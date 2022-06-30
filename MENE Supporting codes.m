% Supporting codes for the MENE survey

% To test and assign Constituencies from MENE into
% standarised named Constituencies
for x = 1:468370
    if ceil (x/100) == x/100
        clc
        fprintf('Percentage of cells done: %2.3f \n',((x/468370)*100));
    end
   num = strfind(PollingCons,char(Constituants(x)));
   counter = 0;
    for i = 1:631 
        if cell2mat(num(i)) == 1 
            List(x,1) = i;
            counter = 1;
        end
    end
    if counter == 0
        Test = find(strcmp(WrongList.Name,char(Constituants(x)))==1);
        if isempty(Test)~= 1
            List(x,1) = WrongList.Number(Test);
        else
            List(x,1) = -888;
        end
    end
    clear counter i num Test
end

ListerNineNine = [];
ErrorList = find(List==-999);
for x = 1:length(ErrorList);
    ListerNineNine(x,1) = ErrorList(x);
end
clear ErrorList
ListerFullyUnknown = [];
ErrorList = find(List==-888);
for x = 1:length(ErrorList);
    ListerFullyUnknown(x,1) = Constituants(ErrorList(x));
end

% Calculate the mean weekly visiting frequency per constituency and number
% of respondends
AverCons = dataset(NaN(631,1),'VarNames','MeanWeekly');
AverCons.Count = NaN(631,1);
TestCount = [];
for i = 1:631
    display(i)
TestCount  = find(Respondents.PolConsNum ==PollingNum(i));
if isempty(TestCount)~=1
AverCons.MeanWeekly(i,1) = nanmean(Respondents.Weekly(TestCount));
else
    AverCons.MeanWeekly(i,1) = NaN;
end
AverCons.Count(i,1) = length(TestCount);
clear test TestCount
end

% To test and assign Constituencies from MENE into
% standarised named Regions
for x = 1:468370
    if ceil (x/100) == x/100
        clc
        fprintf('Percentage of cells done: %2.3f \n',((x/468370)*100));
    end
   num = strfind(CermonialMeme.Name,char(Respondents.Cermonial(x)));
   counter = 0;
    for i = 1:49 
        if cell2mat(num(i)) == 1 
            ListCermonial(x,1) = CermonialMeme.Code(i);
            counter = 1;
        end
    end
    if counter == 0 
        ListCermonial(x,1) = -888;
    end
    clear counter i num Test
end
ListerFullyUnknownCermonial = [];
ErrorListCermonial = find(List==-888);
for x = 1:length(ErrorListCermonial);
    ListerFullyUnknownCermonial(x,1) = Respondents.ObjectID(x));
end

% Calculate the mean weekly visiting frequency per Region and number
% of respondends
AverNuts2 = dataset(NaN(33,1),'VarNames','MeanWeekly');
for i = 1:33
    display(i)
test = strfind(Respondents.CermonialCode,char(NUTS2.Nuts2Code(i)));
count = 1;
TestCount = [];
for j = 1:468370
   tmp = cell2mat(test(j));
   if isempty(tmp)~=1
       TestCount(count) = j;
       count = count +1;
   end
   clear tmp
end
if isempty(TestCount)~=1
AverNuts2.MeanWeekly(i,1) = nanmean(Respondents.Weekly(TestCount));
else
    AverNuts2.MeanWeekly(i,1) = -9999;
end
AverNuts2.Count(i,1) = count-1;
clear test TestCount
end