% spkStats(path,chanmap0ind)
% 
% Get mean spike waveform and calculate various spike statistics. 
% Requires dat files and Kilosort/PhyTemplateGUI sorting results.
% 
% INPUT
%   basepath:    path to your folder containing dat files etc.
%   chanmap0ind: channel map matrix with 0 index. 
%                spike waveform will be plotted according to this map. 
%                Below is an example map for the Buzsaki32 probe.
%                    chanmap0ind = [
%                       16 21 8 13;
%                       19 22 11 14;
%                       18 25 0 15;
%                       17 30 7 12;
%                       20 27 9 2;
%                       28 23 5 10;
%                       24 29 1 4;
%                       26 31 3 6];
% 
% OUTPUT
%   spkStats.mat
%       stats(clu#).
%           id:          original cluster id in Phy
%           maxSh:       maximum-amplitude shank
%           maxCh:       maximum-amplitude channel in the shank (from top to bottom)
%           waveMean:    mean spike waveform of all channels in the shank
%           waveStd:     standard deviation of waveforms
%           waveSem:     standard error of the mean of waveforms
%           troughAmp:   trough amplitude (uV)
%           rise2trough: rise to trough time (ms)
%           trough2peak: trough to peak time (ms)
%           FWHM:        full width at half maximum (ms)
%           isoDist:     isolation distance (well isolated, >20)
%           isiIndex:    ISI index (well isolated, <0.2)
%           spkNum:      number of spikes
%           meanRate:    mean firing rate (Hz)
%       t: time vecor for plotting spike waveforms (ms)
% 
% References
%   Isolation distance: Schmitzer-Torbert et al., Neuroscience, 2005
%   ISI index:          Fee et al, J Neurosci Methods, 1996
% 
% Takuma Kitanishi, OCU, 2017

function spkStats(basepath,chanMap0indMatrix)

% % % test path 
% basepath = 'G:\tk0056\171015PhySorted'
% load('D:\data\rec\config\Kilosort\Buzsaki256.mat','chanMap0indMatrix')


% setting
fprintf('%s%s\n','Gathering spike waveforms: ',basepath)

dat2uV  = 5*10^6/400/(2^16/2);       % dat to micro volt conversion
tmpfile = fullfile(basepath,'spkStatsTmp.mat');

n1 = size(chanMap0indMatrix,1);     % number of channels in a shank
n2 = size(chanMap0indMatrix,2);     % number of shanks
chanmap = chanMap0indMatrix(:)+1;
nmax = length(chanmap(:));

% load spike sample (ss) and clusters (clu)
[ts,clu,cluOri,Fs]=loadSpk(basepath);
% spike samples
ss = round(ts*Fs);

% find a dat file
datfile = dir(fullfile(basepath,'*.dat'));
if length(datfile)<1; error('No dat file detected!'); end
if length(datfile)>1; error('Two or more dat file detected!'); end
bytes   = datfile.bytes;
datfile = datfile.name;
sessName = datfile(1:end-4);

% total recorded channel number
params = loadParamsPy(basepath);
NchanTOT = params.n_channels_dat;

% recording length (sec) calculated from the dat file size
sessLength =  bytes / (2*NchanTOT*Fs); 

% setting for fread
NT1 = round((0.6/1000*Fs));  % '0.6/1000' indicates 0.6ms before spikes
NT2 = round((1.4/1000*Fs));  % '1.4/1000' indicates 1.4ms after spikes
NT  = NT1+NT2+1;             % sample number
t   = (-NT1:NT2)/Fs*1000;    % time vector (ms)

% gather spike waveforms
disp('Gathering spike waveforms...')
nclu = max(clu);
stats = struct;
fid = fopen(fullfile(basepath,datfile),'r');

%% determine the largest amp shank by getting sparsely sampled 100 spikes
for ii = 1:nclu
    % spike timing sample of the cluster
    ssTmp = ss(clu==ii);
    n = length(ssTmp);
    
    if n<100
        nspk = 1:n;
    else
        % nspk = [1:50 n-49:n];                     % initial and last 50 spikes
        nspk = 1:floor(n/100):floor(n/100)*99+1;    % sparse sampling
    end
    wavesTmp = single(nan(NT,NchanTOT,length(nspk)));
    
    for jj=1:length(nspk)
        % get recording data around each spike
        fseek(fid, 2*NchanTOT*(ssTmp(nspk(jj))-NT1), 'bof');
        dat = fread(fid, [NchanTOT NT], '*int16');
        dat = single(dat');

        if size(dat,1)==NT
            % baseline correction
            wavesTmp(:,:,jj) = dat-mean([dat(1,:);dat(end,:)]);
        else
            wavesTmp(:,:,jj) = nan(NT,NchanTOT);
        end
    end
    % trim/sort channels
    wavesTmp = wavesTmp(:,chanmap,:);
    
    % channel number with maximum amplitude
%     [~,maxAmpCh] = max( max(nanmean(wavesTmp,3))-min(nanmean(wavesTmp,3)) );
    meanWav = nanmean(wavesTmp,3);
    [~,maxAmpCh] = max(max(meanWav) - meanWav(NT1,:));
    
    % shank & channel(from top to bottom) number with maximum amplitude
    stats(ii).maxSh = ceil(maxAmpCh/n1);
    stats(ii).maxCh = maxAmpCh-(stats(ii).maxSh-1)*n1;
end


%% gather spikes for each shank

if exist(tmpfile,'file')
    disp('Loading previously stored spkStatsTmp.mat...')
    load(tmpfile,'shank');
else
    shank = struct;
    maxSh = [stats.maxSh];
    
    for ii=1:n2
        % cluster number in this shank
        cluList = find(maxSh==ii);
        cluTmp = clu(any(clu==cluList,2));
        ssTmp  = ss(any(clu==cluList,2));
        n = length(ssTmp);
        fprintf('%s%u%s%u\n','Shank: ',ii,' #spikes: ',n);
        
        if n>1
            % channel range
            chrange = (ii-1)*n1+1:ii*n1;

            % get all spikes in the max amplitude shank
            disp('loading spikes...')
            waves = single(nan(NT,n1,n));
            for jj=1:n
                % get waveform around each spike
                fseek(fid, 2*NchanTOT*(ssTmp(jj)-NT1), 'bof');
                dat = fread(fid, [NchanTOT NT], '*int16');
                dat = single(dat');
                
                if size(dat,1) == NT
                    % baseline correction
                    dat = dat-mean([dat(1,:);dat(end,:)]);
                else
                    dat = nan(NT,NchanTOT);
                end
                % trim/sort channels
                dat = dat(:,chanmap);
                waves(:,:,jj) = dat(:,chrange);
            end
            
            % features for isolation distance
            disp('PCA...')
            fet = nan(n,n1*2);
            % power of spike waveform
            e = sum(waves.^2,1);
            fet(:,1:n1) = squeeze(e)';
            % normaized waveforms
            nwaves = waves./repmat(e,NT,1,1);
            % perform pca and take the 1st component
            for jj=1:n1
                [~,score] = pca( squeeze(permute(nwaves(:,jj,:),[2 3 1])) );
                fet(:,n1+jj) = score(:,1);
            end
            
            shank(ii).waves = waves * dat2uV;   % conversion to uV
            shank(ii).clu = cluTmp;
            shank(ii).fet = fet;
        else
            shank(ii).waves = nan;
            shank(ii).clu   = nan;
            shank(ii).fet   = nan;
        end
    end
    
    clear waves nwaves cluTmp fet ssTmp score;
    
%     disp('Saving tmp file...')
%     save(tmpfile,'shank','-v7.3');
end

fclose(fid);


%% mean, std, sem of spike waveform of each cluster
for ii=1:nclu
    maxSh = stats(ii).maxSh;
    waves = shank(maxSh).waves(:,:,shank(maxSh).clu==ii);
    spkNum = size(waves,3);
    
    stats(ii).waveMean  = nanmean(waves,3);
    stats(ii).waveStd   = nanstd(waves,[],3);
    stats(ii).waveSem   = stats(ii).waveStd / sqrt(spkNum);
    
    % plot waveforms for all channels in the shank
    figure(1); clf
    subplot(121)
    
    x = stats(ii).waveMean;
    e = stats(ii).waveStd;
    pad = repmat(-1*[1:n1],NT,1);
    
    x2 = x(:,stats(ii).maxCh);
    e2 = e(:,stats(ii).maxCh);
    pad2 = -stats(ii).maxCh;
    
    sc = 1/100;
    plot([0 0],[0 -n1-1],':'); hold on
    plot(t,x*sc+pad,'k',t,x2*sc+pad2,'r')
    plot(t,(x2+e2)*sc+pad2,'r:',t,(x2-e2)*sc+pad2,'r:'); hold off
    
    xlabel('Time (ms)')
    ylabel('Channels (channel spacing = 100uV)')
    id = cluOri(clu==ii);
    id = id(1);
    title(sprintf('%s%s%u%s%u%s%u%s%u%s',sessName,'_',ii,' (id',id,') (sh',stats(ii).maxSh,' ch',stats(ii).maxCh,')'),'Interpreter','none')


    % calculate waveform params
    tq = (-NT1:0.1:NT2)/Fs*1000;
    xq = interp1(t,x2,tq,'spline');  % 10-fold upsampling
    
    % trough point
    [trough,troughIdx] = min(xq);
    
    % rise point
    d = diff(xq);
    riseIdx = find(d<mean(d(1:Fs/1000))-0.5*std(d),1);
    rise    = xq(riseIdx);
    
    % peak point
    [pks,locs] = findpeaks(xq(troughIdx:end));
    if ~isempty(pks)
        peakIdx = troughIdx + locs(1) -1;
        peak    = pks(1);
    else
        peakIdx = length(tq);
        peak    = xq(peakIdx);
    end
    
    % FWHM
    idxAll = find(xq-trough/2<0);
    FWHM1idx = idxAll(1);
    FWHM2idx = idxAll(end);
    
    % params
    rise2trough = tq(troughIdx)-tq(riseIdx);    % [ms]
    trough2peak = tq(peakIdx)-tq(troughIdx);    % [ms]
    FWHM = tq(FWHM2idx)-tq(FWHM1idx);           % [ms]
    
    % overall mean rate (Hz)
    nspk = size(waves,3);
    meanRate = nspk/sessLength;
    
    % isolation distance 
    inIdx = shank(maxSh).clu==ii;
    md = mahal(shank(maxSh).fet(~inIdx,:),shank(maxSh).fet(inIdx,:));
    nIn = sum(inIdx);
    nOut= sum(~inIdx);
    if nIn>nOut
        isoDist = nan;
    else
        sortedmd = sort(md);
        isoDist = sortedmd(sum(inIdx));
    end
    
    % isi index (Fee et al., 1996, J Neurosci Methods)
    tsTmp = ts(clu==ii);
    isi = diff(tsTmp) *1000;  % (ms)
    isiIndex = (8/2)*( sum(isi<2)/sum(2<=isi & isi<10) );
    
    % plot the largest amplitude waveform
    subplot(222)
    plot(tq,xq,'-'); hold on
    plot(tq(riseIdx),rise,'ro',tq(troughIdx),trough,'ro',tq(peakIdx),peak,'ro')
    plot([tq(FWHM1idx) tq(FWHM2idx)],[trough/2 trough/2],'r')
    hold off
    xlabel('Time (ms)')
    ylabel('Amplitude (uV)')
    title(sprintf('%s%.2f%s%.2f%s%.2f%s%.1f','Rise2tr=',rise2trough,' Tr2pk=',trough2peak,' FWHM=',FWHM,' Rate=',meanRate))

    subplot(224)
    hist(isi(isi<50),50)
    xlabel('ISI (ms)'); 
    ylabel('# spikes')
    title(sprintf('%s%.3f%s%.0f','ISIidx=',isiIndex,' IsoDist=',isoDist))
    
    drawnow;
    
    % save image
    pngImage = sprintf('%s%u%s','spkStats_',ii,'.png');
    epsImage = sprintf('%s%u%s','spkStats_',ii,'.eps');
    print('-dpng',fullfile(basepath,pngImage));
    print('-depsc','-tiff',fullfile(basepath,epsImage));
    
    % register
    stats(ii).id          = id;              % cluster id in Phy
    stats(ii).troughAmp   = -trough;         % trough amplitude (uV)
    stats(ii).rise2trough = rise2trough;     % rise to trough (ms)
    stats(ii).trough2peak = trough2peak;     % trough to peak (ms)
    stats(ii).FWHM        = FWHM;            % full width at half max (ms)
    stats(ii).isoDist     = isoDist;         % isolation distance
    stats(ii).isiIndex    = isiIndex;        % ISI index
    stats(ii).meanRate    = meanRate;        % mean firing rate (Hz)
    stats(ii).spkNum      = spkNum;          % number of spikes
end

figure(2); clf
subplot(421); hist([stats.isoDist],40);      xlabel('Isolation distance'); title(sessName)
subplot(422); hist([stats.isiIndex],40);     xlabel('ISI index')
subplot(423); hist([stats.troughAmp],40);    xlabel('Trough amplitude (uV)')
subplot(424); hist([stats.rise2trough],40);  xlabel('Rise to trough (ms)')
subplot(425); hist([stats.trough2peak],40);  xlabel('Trough to peak (ms)')
subplot(426); hist([stats.FWHM],40);         xlabel('FWHM (ms)')
subplot(427); hist([stats.meanRate],40);     xlabel('Mean Rate (Hz)')
print('-dpng',fullfile(basepath,'spkStats_s1'));

figure(3); clf
subplot(221); scatter([stats.trough2peak],[stats.FWHM]); xlabel('Trough to peak (ms)'); ylabel('FWHM (ms)'); title(sessName)
subplot(222); scatter([stats.rise2trough],[stats.FWHM]); xlabel('Rise to trough (ms)'); ylabel('FWHM (ms)'); title(sessName)
subplot(223); scatter([stats.trough2peak],[stats.meanRate]); xlabel('Trough to peak (ms)'); ylabel('Mean rate (Hz)')
subplot(224); scatter([stats.FWHM],[stats.meanRate]); xlabel('FWHM (ms)'); ylabel('Mean rate (Hz)')
print('-dpng',fullfile(basepath,'spkStats_s2'));

% unit number map
unitNmap = zeros(n1,n2);
for ii=1:nclu
    unitNmap(stats(ii).maxCh,stats(ii).maxSh) = unitNmap(stats(ii).maxCh,stats(ii).maxSh)+1;
end
figure(4); clf
imagesc(unitNmap); colorbar
pbaspect([2400 1600 1]) % probe size
title(sprintf('%s%s%s%u%s','Unit # map: ',sessName,'  (Total unit #: ',nclu,')')); 
xlabel('Shank'); ylabel('Channel')
print('-dpng',fullfile(basepath,'spkStats_s3'));
print('-depsc','-tiff',fullfile(basepath,'spkStats_s3'));


% save
save(fullfile(basepath,'spkStats.mat'),'stats','t','sessName','chanMap0indMatrix')

%%
% [ts,clu,cluOri,Fs]=loadSpk(basepath,sess)
% 
% Load Kilosort/PhyTemplateGUI clustering. Requires following 3 files: 
%   params.py
%   spike_times.npy
%   spike_clusters.npy
%   cluster_groups.csv
% 
% INPUT
%   basepath:   the folder containing Kilosort clustering results (npy files)
%   sess:       (optional) session number in integer (1,2,3...) or session
%               name (tk0056-171015-06 etc). 
%   
% OUTPUT
%   ts:     spike timing of good (in Phy) clusters [sec]. If session is
%           specified, ts is offset.
%   clu:    sorted good (in Phy) cluster number (1, 2, ...)
%   cluOri: original cluster number generated by Kilosort/PhyTemplateGUI
%   Fs:     sampling frequency [Hz]
% 
% Takuma Kitanishi, OCU, 2017

function [ts,clu,cluOri,Fs]=loadSpk(basepath,varargin)

% No session specified
if nargin==1
    sessFlag = false;
    disp('No session specified. Load the entire period of spike timing.')
end
% Session specified
if nargin==2
    sessFlag = true;
    sess = varargin{1};
    
    % load info file
    [sessName,sessLength] = loadInfo(basepath);

    % find session number
    if isnumeric(sess)
        sessNum = sess;
    elseif ischar(sess)
        sessNum = find(strcmp(sessName,sess));
    else
        error('Unknown input format!')
    end
    
    % session start/end time(sess)
    if isempty(sessNum)
        error('Cannot find the valid session!')
    end
    if sessNum==1
        tsStart = 0;
        tsEnd = sessLength(1);
    else
        cum = cumsum(sessLength);
        tsStart = cum(sessNum-1);
        tsEnd = cum(sessNum);
    end
    fprintf('%s%s%s%.1f%s%.1f%s\n','Session: ',sessName{sessNum},' (',tsStart,'-',tsEnd,' sec)')
end
if nargin>2
    error('Too many inputs!')
end

% check whether previously loaded mat file exists
matFlag = false;
filename = fullfile(basepath,'loadSpk.mat');
if exist(filename,'file')
    % check timestamp of 
    load(filename,'datenum1','datenum2')
    
    info1 = dir(fullfile(basepath,'spike_clusters.npy'));
    info2 = dir(fullfile(basepath,'cluster_groups.csv'));
    if datenum1 == info1.datenum && datenum2 == info2.datenum
        matFlag = true;
    end
end

% load spikes
if matFlag
    disp('Loading spikes from previously saved loadSpk.mat file...')
    load(filename,'ts','clu','cluOri','Fs');
else
    disp('Loading spikes from npy/csv files...')
    
    % get sampling rate from params.py
    params = loadParamsPy(basepath);
    Fs = params.sample_rate;
    
    % load spike timing (unit in samples)
    ssTmp  = npy2mat(fullfile(basepath,'spike_times.npy'));
    % load cluster number
    cluTmp = npy2mat(fullfile(basepath,'spike_clusters.npy'))';
    
    % load clustering result
    % 0=noise, 1=MUA, 2=Good, 3=unsorted
    fid = fopen(fullfile(basepath,'cluster_groups.csv'),'r');
    csv = textscan(fid,'%s %s');
    fclose(fid);
    
    
    id = cellfun(@str2num,csv{1,1}(2:end));
    groupTmp = csv{1,2}(2:end);
    
    n = length(groupTmp);
    group = nan(n,1);
    for ii=1:n
        switch groupTmp{ii}
            case 'noise'
                group(ii) = 0;
            case 'mua'
                group(ii) = 1;
            case 'good'
                group(ii) = 2;
            case 'unsorted'
                group(ii) = 3;
        end
    end
    
    % take only 'good' clusters
    
    % good cluster id
    id_good = id(group==2);
    
    % spike index of 'good' clusters
    idx = logical(sum(cluTmp(:) == id_good',2));
    
    ts  = ssTmp(idx)/Fs;
    cluOri = cluTmp(idx);
    
    % sort cluster number
    [~,~,clu]=unique(cluOri);
    
    % date modified
    info1 = dir(fullfile(basepath,'spike_clusters.npy'));
    datenum1 = info1.datenum;
    info2 = dir(fullfile(basepath,'cluster_groups.csv'));
    datenum2 = info2.datenum;
    
    % save
    save(fullfile(basepath,'loadSpk.mat'),'ts','clu','cluOri','Fs','datenum1','datenum2')
end

% select ts in the session
if sessFlag
    idx = tsStart<ts & ts<tsEnd;
    ts = ts(idx)-tsStart;
    clu = clu(idx);
    cluOri = cluOri(idx);
end


% read other parameters
% amp            = npy2mat(fullfile(path,'amplitudes.npy'));
% chanMap        = npy2mat(fullfile(path,'channel_map.npy'));
% chanPos        = npy2mat(fullfile(path,'channel_positions.npy'));
% pcFetInd       = npy2mat(fullfile(path,'pc_feature_ind.npy'));
% pcFet          = npy2mat(fullfile(path,'pc_features.npy'));
% simTemplates   = npy2mat(fullfile(path,'similar_templates.npy'));
% spkTemplates   = npy2mat(fullfile(path,'spike_templates.npy'));
% templateFetInd = npy2mat(fullfile(path,'template_feature_ind.npy'));
% templateFet    = npy2mat(fullfile(path,'template_features.npy'));
% template       = npy2mat(fullfile(path,'templates.npy'));
% templateInd    = npy2mat(fullfile(path,'templates_ind.npy'));
% templateUnw    = npy2mat(fullfile(path,'templates_unw.npy'));
% whitening      = npy2mat(fullfile(path,'whitening_mat.npy'));

%%
% params = loadParamsPy(basepath)
% 
% INPUT
%   basepath:   the folder containing the params.py file
% OUTPUT
%   params:     structure including the following children;
%               dat_path, n_channels_dat, dtype, offset, sample_rate,
%               hp_filtered
% 
% Takuma Kitanishi, OCU, 2017-2018
function params = loadParamsPy(basepath)
fid = fopen(fullfile(basepath,'params.py'),'r');
txt = textscan(fid,'%s%s%s');

for ii=1:6
    switch txt{1}{ii}
        case 'dat_path'
            params.dat_path = txt{3}{ii}(2:end-1);
        case 'n_channels_dat'
            params.n_channels_dat = str2double(txt{3}{ii});
        case 'dtype'
            params.dtype = txt{3}{ii}(2:end-1);
        case 'offset'
            params.offset = str2double(txt{3}{ii});
        case 'sample_rate'
            params.sample_rate = str2double(txt{3}{ii});
        case 'hp_filtered'
            params.hp_filtered = txt{3}{ii};
    end
end
        
fclose(fid);