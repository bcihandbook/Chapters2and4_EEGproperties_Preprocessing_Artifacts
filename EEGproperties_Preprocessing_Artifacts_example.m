clear all;close all;clc;

%% Add bdfmatlab, biosig toolbox for loading data, and FORCe for artifact 
% removal. Edit paths to fit your own locations
addpath(genpath('./bdfmatlab/'));
addpath(genpath('.FORCe/'));
addpath(genpath('./biosig/'));

%% Useful init
load('chanlocs31Elvira.mat');
SamplingRate = 500; % Hz
NChEEG = 31;

%% Butterworth filter
[b, a] = butter(4, 0.1/(SamplingRate/2), 'high');

%% Find all BDF data
DataPath = ''; % Set here the path to the data on your own machine
BDFfiles = dir([DataPath '/*.bdf']);
BDFfiles = {BDFfiles.name};


%% Load data, extract trials and labels
ArtifactTrials = [];
ArtifactLabels = [];
RestTrials = [];
trial = 0;
for f=1:length(BDFfiles)

    header = readbdfheader([DataPath '/' BDFfiles{f}]);
    data = readbdfdata(header);
    data = data';
    
    %% Keep trigger channel
    trigger = data(:,end);
    
    %% Construct EOG/EMG bipoles
    exgbipoles = [data(:,40)-data(:,37)... % EOGh
        (data(:,39)-data(:,36))/2+(data(:,38)-data(:,35))/2 ... % EOGv
        data(:,35)-data(:,36)... % Cheeks bipole
        data(:,33)-data(:,34)... % Mouth bipole
        data(:,[49:54])]; % Monopolar acc/gyro
    %% Keep EEG channels
    ChannelLabels = {header.Channel.Label};
    ChannelLabels = ChannelLabels(1:32);
    ChannelLabels(25) = []; % Remove EOG
    data = data(:,[1:24 26:32]);
    
    %% Apply high-pass above 0.1 Hz
    %data = filtfilt(b,a,data);
    
    %% Get events
    POS = gettrigger(trigger,1);
    TYP = trigger(gettrigger(trigger,1));
    
    
    for event=1:length(TYP)
    
        if(TYP(event) == 783) % If it is rest
            trial = trial + 1; % There is a rest trial before each artifact trial
            RestTrials(trial,:,:) = data(POS(event):POS(event)+SamplingRate*2-1,:);
            RestTrialsEXG(trial,:,:) = exgbipoles(POS(event):POS(event)+SamplingRate*2-1,:);
        end
        
        if(ismember(TYP(event),[261,1077,1078,1079,1081,1082,1092,1093,1094,...
                1099,1100,1111,1113,1114,1116]))
            ArtifactTrials(trial,:,:) = data(POS(event)+SamplingRate*1:POS(event)+SamplingRate*5-1,:);
            ArtifactLabels(trial) = TYP(event);
            ArtifactTrialsEXG(trial,:,:) = exgbipoles(POS(event)+SamplingRate*1:POS(event)+SamplingRate*5-1,:);
        end
    end
    
end

% DC removal the fast way
ArtifactTrials = ArtifactTrials - permute(repmat(squeeze(mean(ArtifactTrials,2)),1,1,size(ArtifactTrials,2)),[1 3 2]);
RestTrials = RestTrials - permute(repmat(squeeze(mean(RestTrials,2)),1,1,size(RestTrials,2)),[1 3 2]);

%% Compute visualizations without any filtering or artifact removal
computeGATopo(ArtifactTrials, ArtifactLabels, RestTrials, ChannelLabels, locs, 'Fz', 0);


%% Apply CAR
for tr=1:size(ArtifactTrials,1)
    ArtifactTrialsCAR(tr,:,:) = squeeze(ArtifactTrials(tr,:,:)) - ...
        repmat(mean(squeeze(ArtifactTrials(tr,:,:)),2),1, NChEEG);
    RestTrialsCAR(tr,:,:) = squeeze(RestTrials(tr,:,:)) - ...
        repmat(mean(squeeze(RestTrials(tr,:,:)),2),1, NChEEG);
end
computeGATopo(ArtifactTrialsCAR, ArtifactLabels, RestTrialsCAR, ChannelLabels, locs, 'Fz', 10);


%% Apply Surface Laplacian
montage = channels2montage(ChannelLabels);
laplacianMatrix = montage2laplacian(montage,'all');
for tr=1:size(ArtifactTrials,1)
    ArtifactTrialsSP(tr,:,:) = squeeze(ArtifactTrials(tr,:,:))*laplacianMatrix;
    RestTrialsSP(tr,:,:) = squeeze(RestTrials(tr,:,:))*laplacianMatrix;
end
computeGATopo(ArtifactTrialsSP, ArtifactLabels, RestTrialsSP, ChannelLabels, locs, 'Fz', 20);


%% Apply FORCe
for tr=1:size(ArtifactTrials,1)
    disp(['FORCe on trial ' num2str(tr)]);
    for win=1:SamplingRate:size(ArtifactTrials,2)
        win
        cleanEEG = FORCe(squeeze(ArtifactTrials(tr,win:win+SamplingRate-1,:))',...
            SamplingRate, locs, 0);
        ArtifactTrialsFORCe(tr,win:win+SamplingRate-1,:) = cleanEEG';
        if(win+SamplingRate-1 <= size(RestTrials,2))
            cleanEEG = FORCe(squeeze(RestTrials(tr,win:win+SamplingRate-1,:))',...
                SamplingRate, locs, 0);
            RestTrialsFORCe(tr,win:win+SamplingRate-1,:) = cleanEEG';
        end
    end
end
% Since FORCe takes long to compute, it is advised you compute it once and
% save the result in a MAT file. Then load this file and comment out the
% computation for subsequent executions of this script
%load('FORCE_Results.mat');
computeGATopo(ArtifactTrialsFORCe, ArtifactLabels, RestTrialsFORCe, ChannelLabels, locs, 'Fz', 30);


%% Apply regression-based artifact removal
% Compute regression coefficients

% Predictor variables (EXG + constant term)
X = [ones(size(ArtifactTrialsEXG,1)*size(ArtifactTrialsEXG,2)+size(RestTrialsEXG,1)*size(RestTrialsEXG,2),1) ...
    [reshape(ArtifactTrialsEXG,[size(ArtifactTrialsEXG,1)*size(ArtifactTrialsEXG,2) size(exgbipoles,2)]) ...
    ; reshape(RestTrialsEXG,[size(RestTrialsEXG,1)*size(RestTrialsEXG,2) size(exgbipoles,2)])]];

% Output variables (EEG)
Y = [reshape(ArtifactTrials,[size(ArtifactTrials,1)*size(ArtifactTrials,2) NChEEG]) ...
    ; reshape(RestTrials,[size(RestTrials,1)*size(RestTrials,2) NChEEG])];

b = [];
for ch=1:NChEEG
    b = [b regress(Y(:,ch), X)];
end
%b = mvregress(X(1:1000,:), Y(1:1000,:)); 
%tic;b = mvregress(X, Y);toc

for tr=1:size(ArtifactTrials,1)
    ArtifactTrialsREGRESS(tr,:,:) = squeeze(ArtifactTrials(tr,:,:)) - ...
        [ones(size(ArtifactTrialsEXG,2),1) squeeze(ArtifactTrialsEXG(tr,:,:))]*b;
    RestTrialsREGRESS(tr,:,:) = squeeze(RestTrials(tr,:,:)) - ...
        [ones(size(RestTrialsEXG,2),1) squeeze(RestTrialsEXG(tr,:,:))]*b;
end
computeGATopo(ArtifactTrialsREGRESS, ArtifactLabels, RestTrialsREGRESS, ChannelLabels, locs, 'Fz', 40);

%% Load Formula E file
SamplingRate = 512;

%Butterworth filter
[b, a] = butter(4, [1 60]/(SamplingRate/2), 'bandpass');

[data, header] = sload([DataPath '/FormulaE.gdf']);
trigger = data(:,end);
data(:,[13 19 32 65 66 67]) = [];
data = filtfilt(b,a,data);
cardata = data - repmat(mean(data,2),1, size(data,2));
Freqs = [1:1:60];
ChannelLabels = header.Label;
ChannelLabels([13 19 32 65 66 67]) = [];
ShowChannel = 'FZ';
ChanInd = find(strcmp(ChannelLabels, ShowChannel));
chandata = data(find(trigger==2):find(trigger==4), ChanInd);
carchandata = cardata(find(trigger==2):find(trigger==4), ChanInd);
[psd, fr] = pwelch(data(:,ChanInd), 512, 256, Freqs, 512);
[carpsd, fr] = pwelch(cardata(:,ChanInd), 512, 256, Freqs, 512);
figure(100);plot(Freqs, psd, 'b');
figure(200);plot(Freqs, carpsd, 'r');