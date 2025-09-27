function [] = computeGATopo(ArtifactTrials, ArtifactLabels, RestTrials,...
    ChannelLabels, locs, ShowChannel, figureoffset)

figure(figureoffset+1);
RNDTrial = randi(size(ArtifactTrials,1));
RNDChannel = randi(32);
plot(squeeze(ArtifactTrials(RNDTrial,:,RNDChannel)),'b');hold on;plot(squeeze(RestTrials(RNDTrial,:,RNDChannel)),'r');

UniqueArtifacts = unique(ArtifactLabels);
for ua = 1:length(UniqueArtifacts)
    %% Grand averages
    figure(figureoffset+2);
    GA(ua,:,:) = squeeze(mean(ArtifactTrials(ArtifactLabels==UniqueArtifacts(ua),:,:),1));
    plot(squeeze(GA(ua,:,find(strcmp(ChannelLabels,ShowChannel)))));hold on;

    %% Topoplots
    figure(figureoffset+3);
    GATopo(ua,:) = mean(squeeze(abs(GA(ua,:,:))));
    subplot(4,4,ua);topoplot(squeeze(GATopo(ua,:)),locs,'maplimits', [0 50]);colorbar;
    xlabel(num2str(UniqueArtifacts(ua)));
end
figure(figureoffset+2)
GARest = squeeze(mean(RestTrials));
plot(squeeze(GARest(:,find(strcmp(ChannelLabels,ShowChannel)))));
title(['Signal Grand Average at ' ShowChannel]);
hold off;
legend(strsplit(num2str([UniqueArtifacts 783])));

figure(figureoffset+3)
GARestTopo = mean(abs(GARest));
subplot(4,4,16);topoplot(GARestTopo,locs,'maplimits', [0 200]);colorbar;
xlabel('Rest');

%% Compute spectra
Freqs = [1:1:60];
ChanInd = find(strcmp(ChannelLabels,ShowChannel));
for tr=1:size(ArtifactTrials,1)
    [ArtifactSpectrum(tr,:)] = pwelch(squeeze(ArtifactTrials(tr,:,ChanInd)), 500, 250, Freqs, 500);
    [RestSpectrum(tr,:)] = pwelch(squeeze(RestTrials(tr,:,ChanInd)), 500, 250, Freqs, 500);
end

for ua = 1:length(UniqueArtifacts)
    ArtifactSpectrumGA(ua,:) = mean(ArtifactSpectrum(ArtifactLabels==UniqueArtifacts(ua),:),1);
    figure(figureoffset+4);
    plot(squeeze(ArtifactSpectrumGA(ua,:)));hold on;
end
RestSpectrumGA = mean(RestSpectrum,1);
plot(squeeze(RestSpectrumGA));hold on;
title(['Spectrum Grand Average at ' ShowChannel]);
hold off;
legend(strsplit(num2str([UniqueArtifacts 783])));
