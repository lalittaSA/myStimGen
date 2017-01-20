function [td,s] = stimGen_dynamic_HL(loFreq,hiFreq,toneDur,toneSOA,freqType,cohLevel,breakType,breakTime)
% stimGen_static_HL generates a sequence of tone bursts of 2 frequencies (loFreq & hiFreq)
% each tone lasts 'toneDur' ms and is followed by an inter-tone interval of 'toneSOA'
% the sequence of loFreq and hiFreq tones is determined by 'freqType'
% there can be different levels of coherence 'cohLevel'
% there can be noisy / silence interuption depending on 'breakType'

% INPUT:
% loFreq    - lowest frequency in Hz
% hiFreq    - highest frequency in Hz
% toneDur   - tone duration in ms
% toneSOA   - stimulus onset asynchrony = inter-tone interval in ms
% freqType  - frequency types: LLL | LLH | LHL | HLL | HHL | HLH | LHH | HHH
% cohLevel   - coherence level for high-frequency tone (1-cohLevel - for low-frequency tone)
% breakType - information break (N - noise & S - silence & C - conflict): 

% OUTPUT:
% td - time bins
% s - auditory stream

%% some default variable for function testing
% loFreq = 625; %hz      312.5 |  625 | 1250 | 2500 |  5000
% hiFreq = 1250; %hz     625   | 1250 | 2500 | 5000 | 10000
% toneDur = 40; %ms
% toneSOA = 10; %ms
% trialDur = 3000;
% freqType = 'HHH';
% cohLevel = 0.9;
% breakType = 'N'; %'N'|'S'|'C'
% breakTime = 1000:2000;

%% generate sequence of frequencies
% 3 sections of frequencies in one trial -> freqSec
% 30 bins for frequency dynamic

freqSec = 3;
nBin = 30;
subBin = nBin/freqSec;
%% gradual change
binCenters = subBin/2:subBin:nBin;

hiCoh = cohLevel;
loCoh = 1-cohLevel;
noiseCoh = 0.5;

switch freqType
    case 'LLL', coh = loCoh*ones(1,nBin);
        
    case 'LLH'
        rampCoh = loCoh:(hiCoh-loCoh)/10:hiCoh;
        coh = [loCoh*ones(1,binCenters(2)-1) rampCoh hiCoh*ones(1,binCenters(1))];
    case 'LHH'
        rampCoh = loCoh:(hiCoh-loCoh)/10:hiCoh;
        coh = [loCoh*ones(1,binCenters(1)) rampCoh hiCoh*ones(1,binCenters(2)-1)];
    case 'HLL'
        rampCoh = loCoh:(hiCoh-loCoh)/10:hiCoh;
        coh = [hiCoh*ones(1,binCenters(1)) fliplr(rampCoh) loCoh*ones(1,binCenters(2)-1)];
    case 'HHL'
        rampCoh = loCoh:(hiCoh-loCoh)/10:hiCoh;
        coh = [hiCoh*ones(1,binCenters(2)-1) fliplr(rampCoh) loCoh*ones(1,binCenters(1))];
    case 'LHL'
        rampCoh = loCoh:(hiCoh-loCoh)/10:hiCoh;
        coh = [loCoh*ones(1,binCenters(1)-1) rampCoh fliplr(rampCoh) loCoh*ones(1,binCenters(1)-1)];
    case 'HLH'
        rampCoh = loCoh:(hiCoh-loCoh)/10:hiCoh;
        coh = [hiCoh*ones(1,binCenters(1)-1) fliplr(rampCoh) rampCoh hiCoh*ones(1,binCenters(1)-1)];
    
    case 'HHH', coh = hiCoh*ones(1,nBin);
        
    case 'NNN', coh = noiseCoh*ones(1,nBin);   
       
    case 'LNN'
        rampCoh = loCoh:(noiseCoh-loCoh)/10:noiseCoh;
        coh = [loCoh*ones(1,binCenters(1)) rampCoh noiseCoh*ones(1,binCenters(2)-1)];
    case 'NLN'
        rampCoh = loCoh:(noiseCoh-loCoh)/10:noiseCoh;
        coh = [noiseCoh*ones(1,binCenters(1)-1) fliplr(rampCoh) rampCoh noiseCoh*ones(1,binCenters(1)-1)];
    case 'NNL'
        rampCoh = loCoh:(noiseCoh-loCoh)/10:noiseCoh;
        coh = [noiseCoh*ones(1,binCenters(2)-1) fliplr(rampCoh) loCoh*ones(1,binCenters(1))];
        
    case 'LLN'
        rampCoh = loCoh:(noiseCoh-loCoh)/10:noiseCoh;
        coh = [loCoh*ones(1,binCenters(2)-1) rampCoh noiseCoh*ones(1,binCenters(1))];
    case 'LNL'
        rampCoh = loCoh:(noiseCoh-loCoh)/10:noiseCoh;
        coh = [loCoh*ones(1,binCenters(1)-1) rampCoh fliplr(rampCoh) loCoh*ones(1,binCenters(1)-1)];
    case 'NLL'
        rampCoh = loCoh:(noiseCoh-loCoh)/10:noiseCoh;
        coh = [noiseCoh*ones(1,binCenters(1)) fliplr(rampCoh) loCoh*ones(1,binCenters(2)-1)];

    case 'HNN'
        rampCoh = noiseCoh:(hiCoh-noiseCoh)/10:hiCoh;
        coh = [hiCoh*ones(1,binCenters(1)) fliplr(rampCoh) noiseCoh*ones(1,binCenters(2)-1)];
    case 'NHN'
        rampCoh = noiseCoh:(hiCoh-noiseCoh)/10:hiCoh;
        coh = [noiseCoh*ones(1,binCenters(1)-1) rampCoh fliplr(rampCoh) noiseCoh*ones(1,binCenters(1)-1)];
    case 'NNH'
        rampCoh = noiseCoh:(hiCoh-noiseCoh)/10:hiCoh;
        coh = [noiseCoh*ones(1,binCenters(2)-1) rampCoh hiCoh*ones(1,binCenters(1))];
        
    case 'HHN'
        rampCoh = noiseCoh:(hiCoh-noiseCoh)/10:hiCoh;
        coh = [hiCoh*ones(1,binCenters(2)-1) fliplr(rampCoh) noiseCoh*ones(1,binCenters(1))];
    case 'HNH'
        rampCoh = noiseCoh:(hiCoh-noiseCoh)/10:hiCoh;
        coh = [hiCoh*ones(1,binCenters(1)-1) fliplr(rampCoh) rampCoh hiCoh*ones(1,binCenters(1)-1)];
    case 'NHH'
        rampCoh = noiseCoh:(hiCoh-noiseCoh)/10:hiCoh;
        coh = [noiseCoh*ones(1,binCenters(1)) rampCoh hiCoh*ones(1,binCenters(2)-1)];
        
    case 'HLN'
        rampCoh1 = loCoh:(hiCoh-loCoh)/10:hiCoh;
        rampCoh2 = loCoh:(noiseCoh-loCoh)/10:noiseCoh;
        coh = [hiCoh*ones(1,binCenters(1)-1) fliplr(rampCoh1) rampCoh2 noiseCoh*ones(1,binCenters(1)-1)];
    case 'LHN'
        rampCoh1 = loCoh:(hiCoh-loCoh)/10:hiCoh;
        rampCoh2 = noiseCoh:(hiCoh-noiseCoh)/10:hiCoh;
        coh = [loCoh*ones(1,binCenters(1)-1) rampCoh1 fliplr(rampCoh2) noiseCoh*ones(1,binCenters(1)-1)];
    case 'HNL'
        rampCoh1 = noiseCoh:(hiCoh-noiseCoh)/10:hiCoh;
        rampCoh2 = loCoh:(noiseCoh-loCoh)/10:noiseCoh;
        coh = [hiCoh*ones(1,binCenters(1)-1) fliplr(rampCoh1) fliplr(rampCoh2) loCoh*ones(1,binCenters(1)-1)];
    case 'LNH'
        rampCoh1 = loCoh:(noiseCoh-loCoh)/10:noiseCoh;
        rampCoh2 = noiseCoh:(hiCoh-noiseCoh)/10:hiCoh;
        coh = [loCoh*ones(1,binCenters(1)-1) rampCoh1 rampCoh2 hiCoh*ones(1,binCenters(1)-1)];
    case 'NLH'
        rampCoh1 = loCoh:(noiseCoh-loCoh)/10:noiseCoh;
        rampCoh2 = loCoh:(hiCoh-loCoh)/10:hiCoh;
        coh = [noiseCoh*ones(1,binCenters(1)-1) fliplr(rampCoh1) rampCoh2 hiCoh*ones(1,binCenters(1)-1)];
    case 'NHL'
        rampCoh1 = noiseCoh:(hiCoh-noiseCoh)/10:hiCoh;
        rampCoh2 = loCoh:(hiCoh-loCoh)/10:hiCoh;
        coh = [noiseCoh*ones(1,binCenters(1)-1) rampCoh1 fliplr(rampCoh2) loCoh*ones(1,binCenters(1)-1)];

    otherwise, disp('invalid coherence setting!')
end

%% abrupt change
% switch freqType
%     case 'LLL', coh = zeros(1,nBin);
%     case 'LLH', coh = [zeros(1,subBin*2) ones(1,subBin)];
%     case 'LHH', coh = [zeros(1,subBin) ones(1,subBin*2)];
%     case 'HLL', coh = [ones(1,subBin) zeros(1,subBin*2)];
%     case 'HHL', coh = [ones(1,subBin*2) zeros(1,subBin)];
%     case 'LHL', coh = [zeros(1,subBin) ones(1,subBin) zeros(1,subBin)];
%     case 'HLH', coh = [ones(1,subBin) zeros(1,subBin) ones(1,subBin)];
%     case 'HHH', coh = ones(1,nBin);
%     otherwise, disp('invalid coherence setting!')
% end
%%
numToneinBloc = 2;
toneDurTot = toneSOA+toneDur;
trialDur = nBin * numToneinBloc * toneDurTot;

breakBins = unique(round(breakTime/(toneDurTot*numToneinBloc)));
switch breakType
    case 'none' 
    case 'N',  coh(breakBins) = 0.5;
    case 'S',  coh(breakBins) = nan;
    case 'C'
        if strcmp(cohType,'HHH'), coh(breakBins) = loFreq; end
        if strcmp(cohType,'LLL'), coh(breakBins) = hiFreq; end
    otherwise, disp('invalid information break type!')
end

numCh = 1;
% calc no. of tone presentations
bufferLen = ceil(trialDur/toneDurTot);
numTonePresentations = floor(trialDur/toneDurTot);
frequencyBuffer = ones(bufferLen,numCh)*hiFreq;

stateList = zeros(numTonePresentations,numCh);
for i = 1:numTonePresentations
    stateList(i,:) = randperm(numCh);
end

for i = 1:numCh
    for bb = 1:length(coh)
        ind = (bb-1)*numToneinBloc;
        if ~isnan(coh(bb))
            numToShuffle = round((1-coh(bb))*numToneinBloc);
            stateIdx = randperm(numToneinBloc);
            stateIdx = stateIdx(1,1:numToShuffle);

            stateList(ind+stateIdx,i) = 1;
            
            frequencyBuffer(ind+stateIdx,i) = loFreq; %low
        else
            frequencyBuffer(ind+1:ind+numToneinBloc,i) = nan;
        end
    end
    numLo = sum(frequencyBuffer(:,i) == loFreq);
    numHi = sum(frequencyBuffer(:,i) == hiFreq);
end

%% generate sequence of tones
fs = 348000;
rampLength = 5;
cosRamp = (rampLength/1000)*fs;

multiplier = fs/1000;
s = zeros(1,trialDur*multiplier);
f = nan(1,trialDur*multiplier);
for ii = 1:numTonePresentations
    ind = (ii-1)*(toneSOA+toneDur)*multiplier + 1;
    if ~isnan(frequencyBuffer(ii))
        tmax = toneDur/1000;                               % Time Duration Of Signal (sec)
        t = linspace(0, tmax, tmax*fs);                          % Time Vector
        freq = frequencyBuffer(ii);                                 % Original Frequency
        tone = sin(2*pi*freq*t);% Original Signal %TODO: multiply voltage - from calib
        tone = pa_ramp(tone, cosRamp, fs);
        s(ind:ind+toneDur*multiplier-1) = tone';
        f(ind:ind+toneDur*multiplier-1) = freq;
    end
end

td = 0:1/fs:(length(s)-1)/fs;

%% for testing - listen
% sound(s, fs) 

% record
% filename = 'sample_stim2.wav';
% audiowrite(filename,s,fs)
%% for testing - plot
% temp_fig_path = '/dataAnalysis/git_public/Penn_auditoryDecision/stimuli/';
% fig_name = 'sample_stim2_toneplot';
% h = figure('Name',fig_name,'Position',get(0,'ScreenSize'));
% subplot(2,1,1)
% tone = sin(2*pi*freq*t);
% plot(tone)
% tone = pa_ramp(tone, cosRamp, fs);
% subplot(2,1,2)
% plot(tone)
% saveas(h,[temp_fig_path fig_name '.fig']);
% save2pdf([temp_fig_path fig_name '.pdf'],h);
% clf
% 
% % 
% % plot(td,s)
% % xlabel('Time (ms)')
% % ylabel('Amplitude')
% % title('Signal')
% 
% temp_fig_path = '/dataAnalysis/git_public/Penn_auditoryDecision/stimuli/';
% fig_name = 'sample_stim2_plot';
% h = figure('Name',fig_name,'Position',get(0,'ScreenSize'));
% subplot(2,1,1)
% plot(td,f/1000)
% text(td(end),(hiFreq+100)/1000,['numHi = ' num2str(numHi)],'HorizontalAlignment','right');
% text(td(end),(loFreq-100)/1000,['numLo = ' num2str(numLo)],'HorizontalAlignment','right');
% xlabel('Time (secs)')
% ylabel('Frequency(Hz)')
% ylim([0 3])
% 
% subplot(2,1,2)
% segmentLength = round(numel(s)/5);
% % window = hamming(512);
% % noverlap = 256;
% % nfft = 512;
% % [S,F,T,P] = spectrogram(s,window,noverlap,nfft,fs,'yaxis');
% % surf(T,F,10*log10(P),'edgecolor','none'); axis tight; view(0,90);
% % colormap(hot)
% % set(gca,'clim',[-80 -30])
% % ylim([0 3])
% % 
% % periodogram(s,[],[],fs)
% 
% % spectrogram(s,segmentLength,[],[],fs,'yaxis')
% spectrogram(s,round(segmentLength/5),round(80/100*segmentLength/5),[],fs,'yaxis')
% % spectrogram(s,round(numel(s)/25),[],round(numel(s)/25),fs,'yaxis')
% ylim([0 3])
% saveas(h,[temp_fig_path fig_name '.fig']);
% save2pdf([temp_fig_path fig_name '.pdf'],h);
% caf