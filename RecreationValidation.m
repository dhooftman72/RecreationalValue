function RecreationValidation
warning off

for Validatorset = 1:1:2
    clearvars -except Validatorset
    [Validators,Pressure,Transfer,DistoLondon] = DefineParameters(Validatorset);
    
    NrValidators = length( Transfer.NamelistValidators);
    NrDependents = length( Transfer.NamelistPressure);
    Transfer.ToRemove = []; % since I have not found a dataset alternative for 'exist'
    Transfer.ValidLength = []; % idem
    Dummy = NaN(NrDependents,NrValidators);
    Ranking.RHO =dataset(Dummy(:,1),'ObsNames',Transfer.NamelistPressure,'VarNames',Transfer.NamelistValidators(:,1));
    for i = 2:NrValidators
        Ranking.RHO.(genvarname(char(Transfer.NamelistValidators(i)))) = Dummy(:,i);
    end
    Ranking.PVAL = Ranking.RHO;
    Mean_double_deviation = Ranking.RHO;
    Variables.Fvalue = Ranking.RHO;
    Variables.Pvalue = Ranking.RHO;
    Variables.Rico = Ranking.RHO;
    Variables.Rsquared = Ranking.RHO;
    DistancetoLondon.FValue = dataset((Dummy(1,1)'),'VarNames',Transfer.NamelistPressure(1));
    for i = 2:NrDependents
        DistancetoLondon.FValue.(genvarname(char(Transfer.NamelistPressure(i)))) = Dummy(1,i);
    end
    DistancetoLondon.PValue = DistancetoLondon.FValue;
    DistancetoLondon.Rico = DistancetoLondon.FValue;
    DistancetoLondon.Rsquared = DistancetoLondon.FValue;
    clear i Dummy
    
    for Vali = 1:NrValidators
        display('  ')
        display(Transfer.NamelistValidators(Vali))
        for depen = 1:NrDependents
            Testset = dataset(Validators(:,Vali),'Varnames','Validator');
            Testset.Pressure = Pressure(:,depen);
            Testset.DistVar = DistoLondon;
            %Make sure No data is No data
            Testset.Validator(Testset.Validator <0) = NaN;
            Testset.Pressure(Testset.Pressure <0) = NaN;
            Testset.DistVar(Testset.DistVar <0) = NaN;
            % Clean non existend values
            clear x_range y_range
            [Transfer] = CleanOutNaN(Testset,Transfer);
            Testset(Transfer.ToRemove,:) = [];
            % Normalise for Deviation, Note before bootstrap
            Validator_range = WinsorFunction(Testset.Validator,2.5);
            Dependent_range = WinsorFunction(Testset.Pressure,2.5);
            DisLondon_range = WinsorFunction(Testset.DistVar,2.5);
            
            %% For under and over use
            if Vali >= 8
                % Ranking difference
                RankVali = tiedrank(Testset.Validator);
                RankY = tiedrank(Testset.Pressure);
                Ranksingle(:,1) = (((RankVali-RankY)))/(length(Testset.Validator)./3);
                % 1) the expected difference under a random is 1/3 times length
                % (I have tested that with 1,000,000 random draws)
                % 2) If the ranking difference increases (so less accuracy) the value
                % of Dsingle increases, so ranking error is good (low) to
                % bad(high)
                Overuse.Ranking.(genvarname(char(Transfer.NamelistValidators(Vali))))(Transfer.OriginalOrder,depen) =...
                    (Ranksingle(:,1)./(max(Ranksingle(:,1))- min(Ranksingle(:,1))));
                Overuse.Ranking.(genvarname(char(Transfer.NamelistValidators(Vali))))(Transfer.ToRemove,depen) = -9999;
                clear RankVali RankY Ranksingle
                
                % Correct for mean value and std differences of normalised
                % values: so they have the same mean and std
                [ValidatorAdapted,DependentAdapted] = SynchroDistri(Validator_range,Dependent_range);
                %Calculate proportional overuse on mean & std corrected
                % difference (so they have the same mean and std)
                OveruseProportional(Transfer.OriginalOrder,1) = ValidatorAdapted - DependentAdapted;  %#ok<*AGROW>
                OveruseProportional(Transfer.ToRemove) = NaN;
                OveruseProportional = OveruseProportional - nanmean(OveruseProportional);
                OveruseProportional(OveruseProportional<-1) = -1;
                OveruseProportional(OveruseProportional>1) = 1;
                OveruseProportional(isnan(OveruseProportional)==1) = -9999;
                Overuse.Proportional.(genvarname(char(Transfer.NamelistValidators(Vali))))(:,depen) = OveruseProportional;
            end
            %% Bootstrap
            % Initiate bootstrap files
            tmpRanking.RHO = NaN(Transfer.bootmax,1);
            tmpRanking.PVAL = NaN(Transfer.bootmax,1);
            tmpMean_double_deviation  = NaN(Transfer.bootmax,1);
            tmpVariables.Fvalue = NaN(Transfer.bootmax,1);
            tmpVariables.Pvalue = NaN(Transfer.bootmax,1);
            tmpVariables.Rico = NaN(Transfer.bootmax,1);
            tmpVariables.Rsquared = NaN(Transfer.bootmax,1);
            tmpDistancetoLondon.PValue = NaN(Transfer.bootmax,1);
            tmpDistancetoLondon.FValue = NaN(Transfer.bootmax,1);
            tmpDistancetoLondon.Rico = NaN(Transfer.bootmax,1);
            tmpDistancetoLondon.Rsquared = NaN(Transfer.bootmax,1);
            % Start bootstrap
            if Validatorset == 1
                display('  ')
            end
            for boot = 1:1:Transfer.bootmax
                if boot == 1
                    fprintf('Doing dependent: %s \n', char(Transfer.NamelistPressure(depen)))
                elseif boot == 2 || (ceil(boot/25)== (boot/25))
                    fprintf('Doing bootstrap nr %i from total %i for %s \n', boot,Transfer.bootmax,char(Transfer.NamelistPressure(depen)))
                end
                % Perform bootstrap
                tmper = randperm(Transfer.ValidLength);
                tmperList = sort(tmper(1:Transfer.DFtarget));
                x_range = Validator_range(tmperList);
                y_range =  Dependent_range(tmperList);
                Dis_range = DisLondon_range(tmperList);
                %Spearman rank
                [tmpRanking.RHO(boot)] = corr(x_range,y_range,'type','Spearman');
                [~,tmpRanking.PVAL(boot)] = corr(x_range,y_range,'type','Spearman');
                % Inversed Deviance
                Datapoints = size(x_range,1);
                Deviation_point= abs(y_range-x_range); %% Accuracy per point
                tmpMean_double_deviation(boot) = 1- (nansum(Deviation_point)/Datapoints); %%Accuracy overall
                %Mean_double_deviation(depen,Vali) = (round(Mean_double.*1000))./1000;
                clear  Deviation_point  Mean_double Datapoints
                
                % Anova on normalised data Distance to London R2.
                [~,outs,statsDis] = anovan(y_range,{Dis_range},'sstype',1,'model',[1],...
                    'continuous', [1], 'display', 'off','varnames', {'DistancetoLondon'});
                Points_Corrected =  y_range - statsDis.resid;
                tmpDistancetoLondon.FValue(boot) = cell2mat(outs(2,6));
                tmpDistancetoLondon.PValue(boot) = cell2mat(outs(2,7));
                tmpDistancetoLondon.Rico(boot) = statsDis.coeffs(2,1);
                tmpDistancetoLondon.Rsquared(boot) = cell2mat(RsquaredFunc(Points_Corrected,y_range));
                clear outs
                % Anova on normalised data for Variabele
                [~,outs,stats] = anovan(statsDis.resid,{x_range},'sstype',1,...
                    'model',[1],'continuous', [1], 'display', 'off',...
                    'varnames', {char(Transfer.NamelistValidators(Vali))});
                Points_Corrected2 =  statsDis.resid - stats.resid;
                tmpVariables.Fvalue(boot) = cell2mat(outs(2,6));
                tmpVariables.Pvalue(boot) = cell2mat(outs(2,7));
                tmpVariables.Rico(boot) = stats.coeffs(2,1);
                tmpVariables.Rsquared(boot) = cell2mat(RsquaredFunc(Points_Corrected2,statsDis.resid));
                clear Testset x_range y_range Dis_range DistVar outs Points_Corrected* stats*
            end
            %% collate Bootstraps
            Ranking.RHO.(genvarname(char(Transfer.NamelistValidators(Vali))))(depen,1) = nanmean(tmpRanking.RHO);
            Ranking.PVAL.(genvarname(char(Transfer.NamelistValidators(Vali))))(depen,1) = nanmean(tmpRanking.PVAL);
            Mean_double_deviation.(genvarname(char(Transfer.NamelistValidators(Vali))))(depen,1) = nanmean(tmpMean_double_deviation);
            Variables.Fvalue.(genvarname(char(Transfer.NamelistValidators(Vali))))(depen) = nanmean(tmpVariables.Fvalue);
            Variables.Pvalue.(genvarname(char(Transfer.NamelistValidators(Vali))))(depen) = nanmean(tmpVariables.Pvalue);
            Variables.Rico.(genvarname(char(Transfer.NamelistValidators(Vali))))(depen) = nanmean(tmpVariables.Rico);
            Variables.Rsquared.(genvarname(char(Transfer.NamelistValidators(Vali))))(depen) = nanmean(tmpVariables.Rsquared);
            DistancetoLondon.PValue.(genvarname(char(Transfer.NamelistPressure(depen)))) = nanmean(tmpDistancetoLondon.PValue);
            DistancetoLondon.FValue.(genvarname(char(Transfer.NamelistPressure(depen)))) = nanmean(tmpDistancetoLondon.FValue);
            DistancetoLondon.Rico.(genvarname(char(Transfer.NamelistPressure(depen)))) = nanmean(tmpDistancetoLondon.Rico);
            DistancetoLondon.Rsquared.(genvarname(char(Transfer.NamelistPressure(depen)))) = nanmean(tmpDistancetoLondon.Rsquared);
            clear tmp*
        end
    end
    %% Outputs
    InPuts.Validators = Validators; %#ok<*STRNU>
    Inputs.Pressure = Pressure;
    save(Transfer.Outfile,'Ranking','Mean_double_deviation','DistancetoLondon','Variables','Inputs','Overuse');
    clear Ranking Mean_double_deviation DistancetoLondon Variables Inputs Overuse
end
end

function  [OutVar] = WinsorFunction(InVar,percLow)
InVar = reshape(InVar,((size(InVar,1)).*(size(InVar,2))),1);
InVar(InVar<0) = NaN;
prct =  prctile(InVar,percLow);
if prct < 0
    display('zero Values present; correct this first')
    cccc
end
InVar_norm = InVar - prct;
clear InVar
InVar_norm( InVar_norm<0) = 0;
prct(2) = prctile(InVar_norm,(100-percLow));
if exist('InVar_org') == 1; %#ok<EXIST>
    InVar_norm = InVar_org;
end
OutVar = (InVar_norm./prct(2));
clear InVar_norm InVar_org testlist testtmp upboud
OutVar(OutVar>1) = 1;
end

function Var = Normalise(Var,Min,Max)
Var = Var./max(Var);
Var(Var>Max) = Max;
Var(Var <Min) = Min;
end

function Rsquared = RsquaredFunc(ExpecIn,ObservedIn)
ExpecIn = reshape(ExpecIn,[],1);
ObservedIn = reshape(ObservedIn,[],1);
List = find((isnan(ExpecIn)== 1));
List2 = find((isnan(ObservedIn)== 1));
List = [(reshape(List,[],1));(reshape(List2,[],1))];
ExpecIn(List) = [];
ObservedIn(List) = [];
meanY = mean(ObservedIn);
Size = length(ExpecIn);
for t = 1:1:Size
    ssres(t) = ((ExpecIn(t)-ObservedIn(t)).^2);
    sstot(t) =  ((ObservedIn(t)-meanY).^2);
end
Rsquared= {1- ((sum(ssres)) /(sum(sstot)))};
end

function [Transfer] = CleanOutNaN(ArrayIn,Transfer)
if isnumeric(ArrayIn) == 0
    ArrayIn = double(ArrayIn);
else
    cccc
end
Transfer.OriginalOrder = reshape((1:1:length(ArrayIn(:,1))),[],1);
Transfer = rmfield(Transfer, 'ToRemove');
Transfer = rmfield(Transfer, 'ValidLength');
ArrayIn(ArrayIn<0) = NaN;
a1=find((isnan(ArrayIn(:,1))==1));
b=find((isnan(ArrayIn(:,2))==1));
c=find((isnan(ArrayIn(:,3))==1));
all = [a1;b;c];
Transfer.ToRemove(:,1) = unique(all);
Transfer.ValidLength = length(ArrayIn(:,1)) - length(Transfer.ToRemove);
if Transfer.ValidLength < Transfer.DFtarget
    Transfer.DFtarget = Transfer.ValidLength;
end
Transfer.OriginalOrder(Transfer.ToRemove,:) = [];
end

function [Var1Adapted,Var2Adapted] = SynchroDistri(VarIn1,VarIn2)
MeanVar1 = nanmean(VarIn1);
MeanVar2  = nanmean(VarIn2);
StdVar1 = nanstd(VarIn1);
StdVar2 = nanstd(VarIn2);
ZrangeVar1 = (VarIn1 - MeanVar1)./StdVar1;
ZrangeVar2 = (VarIn2 - MeanVar2)./StdVar2;
Meantoshift = (MeanVar1+MeanVar2)./2;
StdtoShift = (StdVar1+StdVar2)./2;
StdRatio.Var1 = StdVar1.*(StdtoShift./StdVar1);
StdRatio.Var2 = StdVar2.* (StdtoShift./StdVar2);
Var1Adapted = Meantoshift + (ZrangeVar1.*StdRatio.Var1);
Var2Adapted = Meantoshift + (ZrangeVar2.*StdRatio.Var2);
Var1Adapted = Normalise(Var1Adapted,0,1);
Var2Adapted = Normalise(Var2Adapted,0,1);
end


function [Validators,Pressure,Transfer,DistoLondon] =...
    DefineParameters(Validatorset)
if Validatorset == 1
    load('ConstituencyData.mat')
    Validators(:,1) = log10(ConstituencyData.Pricem2+1);
    Validators(:,2) = ConstituencyData.PotentialServices./ConstituencyData.Hectares;
    Validators(:,3) = ConstituencyData.RealisedServices./ConstituencyData.Hectares;
    Validators(:,4) = ConstituencyData.MeanAttractiveness;
    Validators(:,5) = ConstituencyData.AreaRatioAttractive;
    Validators(:,6) = ConstituencyData.PathsRation;
    Validators(:,7) = log10(ConstituencyData.DistancetoLondon);
    Validators(:,8) = ConstituencyData.NrPandN;
    Validators(:,9) = ConstituencyData.MemeWeekly;
    Validators(:,10) = ConstituencyData.MemeMonthly;
    Validators(:,11) = ConstituencyData.MemeCount;
    Pressure(:,1) = log10(ConstituencyData.PressureYearlyPaths+1);
    Pressure(:,2) = log10(ConstituencyData.PressureMontlyPaths+1);
    Pressure(:,3) = log10(ConstituencyData.PressureWeeklyPaths+1);
    Pressure(:,4) = log10(ConstituencyData.PressureYearlyAllArea+1);
    Pressure(:,5) = log10(ConstituencyData.PressureMontlyAllArea+1);
    Pressure(:,6) = log10(ConstituencyData.PressureWeeklyAllArea+1);
    Pressure(:,7) = log10((ConstituencyData.PressureYearlyPaths./ConstituencyData.Hectares)+1);
    Pressure(:,8) = log10((ConstituencyData.PressureMontlyPaths./ConstituencyData.Hectares)+1);
    Pressure(:,9) = log10((ConstituencyData.PressureWeeklyPaths./ConstituencyData.Hectares)+1);
    Pressure(:,10) =  ConstituencyData.PathsRation;
    Pressure(:,11) = log10(ConstituencyData.Pricem2+1);
    DistoLondon = log10(ConstituencyData.DistancetoLondon+1);
    
    Transfer.NamelistValidators = [{'Pricem2'},{'PotentialServices'},{'RealisedServices'},...
        {'MeanAttractiveness'},{'AreaRatioAttractive'},{'PathsperArea'},{'DistancetoLondon'},...
        {'NrPandN'},{'MemeWeekly'},{'MemeMonthly'},{'MemeNrPeople'}];
    Transfer.NamelistPressure = [{'PressureYearlyPaths'},{'PressureMontlyPaths'},{'PressureWeeklyPaths'},...
        {'PressureYearlyAllArea'},{'PressureMontlyAllArea'},{'PressureWeeklyAllArea'},...
        {'DensityYearlyPaths'},{'DensityMontlyPaths'},{'DensityWeeklyPaths'},...
        {'PathsperArea'},{'Pricem2'}];
    Transfer.Outfile = 'ConstituenciesResults.mat';
    Transfer.bootmax = 250;
    Transfer.DFtarget = 50;
    
elseif Validatorset == 2
    load('NUTS2RegionData.mat')
    Validators(:,1) = log10(NUTS2RegionData.Pricem2+1);
    Validators(:,2) = NUTS2RegionData.PotentialServices./NUTS2RegionData.Hectares;
    Validators(:,3) = NUTS2RegionData.RealisedServices./NUTS2RegionData.Hectares;
    Validators(:,4) = NUTS2RegionData.MeanAttractiveness;
    Validators(:,5) = NUTS2RegionData.AreaRatioAttractive;
    Validators(:,6) = NUTS2RegionData.PathsRation;
    Validators(:,7) = log10(NUTS2RegionData.DistancetoLondon);
    Validators(:,8) = NUTS2RegionData.NrPandN;
    Validators(:,9) = NUTS2RegionData.MemeWeekly;
    Validators(:,10) = NUTS2RegionData.MemeMonthly;
    Validators(:,11) = NUTS2RegionData.MemeCount;
    Validators(:,12) = NUTS2RegionData.TotalDomestincTourism;
    Pressure(:,1) = log10(NUTS2RegionData.PressureYearlyPaths+1);
    Pressure(:,2) = log10(NUTS2RegionData.PressureMontlyPaths+1);
    Pressure(:,3) = log10(NUTS2RegionData.PressureWeeklyPaths+1);
    Pressure(:,4) = log10(NUTS2RegionData.PressureYearlyAllArea+1);
    Pressure(:,5) = log10(NUTS2RegionData.PressureMontlyAllArea+1);
    Pressure(:,6) = log10(NUTS2RegionData.PressureWeeklyAllArea+1);
    Pressure(:,7) = log10((NUTS2RegionData.PressureYearlyPaths./NUTS2RegionData.Hectares)+1);
    Pressure(:,8) = log10((NUTS2RegionData.PressureMontlyPaths./NUTS2RegionData.Hectares)+1);
    Pressure(:,9) = log10((NUTS2RegionData.PressureWeeklyPaths./NUTS2RegionData.Hectares)+1);
    Pressure(:,10) =  NUTS2RegionData.PathsRation;
    Pressure(:,11) = log10(NUTS2RegionData.Pricem2+1);
    DistoLondon = log10(NUTS2RegionData.DistancetoLondon+1);
    Transfer.NamelistValidators = [{'Pricem2'},{'PotentialServices'},{'RealisedServices'},...
        {'MeanAttractiveness'},{'AreaRatioAttractive'},{'PathsperArea'},{'DistancetoLondon'},...
        {'NrPandN'},{'MemeWeekly'},{'MemeMonthly'},{'MemeNrPeople'},{'TotalDomestincTourism'}];
    Transfer.NamelistPressure = [{'PressureYearlyPaths'},{'PressureMontlyPaths'},{'PressureWeeklyPaths'},...
        {'PressureYearlyAllArea'},{'PressureMontlyAllArea'},{'PressureWeeklyAllArea'},...
        {'DensityYearlyPaths'},{'DensityMontlyPaths'},{'DensityWeeklyPaths'},...
        {'PathsperArea'},{'Pricem2'}];
    Transfer.Outfile = 'NUTS2RegionResults.mat';
    Transfer.bootmax = 1;
    Transfer.DFtarget = 33;
end
end