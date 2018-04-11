% spkStats(basepath,chanMap,options)
%
% Get mean spike waveform and calculate various spike statistics.
% Requires dat files and Kilosort/PhyTemplateGUI sorting results.
%
% INPUT
%   basepath:    path to your folder containing param.py and *.npy files.
%   chanMap: channel map matrix or cell with 0 index.
%            spike waveform will be plotted according to this map.
%
%            Below is an example map for the Buzsaki32 probe on amplepex.
%            chanMap = [
%                       16 21 8 13;
%                       19 22 11 14;
%                       18 25 0 15;
%                       17 30 7 12;
%                       20 27 9 2;
%                       28 23 5 10;
%                       24 29 1 4;
%                       26 31 3 6];
%
%            Below is an example map for the Buzsaki64sp probe on aligned dat.
%            chanMap = {0:9
%                       10:19
%                       20:29
%                       30:39
%                       40:49
%                       50:59
%                       60:63};
%
%  options: pairs of name and value 
%       inventory of option names and default values
%        dat2uV  (5*10^6/400/(2^16/2));      scale functor for conversion from dat to micro volt 
%        nspkMax (100)   number of spikes to determine on which channel the cluster is
%        baselineMethod ('mean') mean or interpolate. Determine how to extact baseline from each spike 
%        excludeNoise (true)     whether exlcude noise from isoDist and Lratio
%        outputFigure (true)     output summary figure for each unit
%        outputOnlyGood (true)   save figures only for good units 
%                                (regardless of setting, calculation will be done for all units) 
%        figSaveDir (basepath)   path to save figures if output figure will be saved
%        matSaveDir (basepath)   path to save final output (spkStats.mat)
%
% OUTPUT
%   spkStats.mat
%       stats(clu#).
%           id:          original cluster id in Phy
%           maxSh:       maximum-amplitude shank
%           maxCh:       maximum-amplitude channel in the shank (from top to bottom)
%           group:       group assigned on Phy (good/mua/noise/unsorted)
%           waveMean:    mean spike waveform of all channels in the shank
%           waveStd:     standard deviation of waveforms
%           waveSem:     standard error of the mean of waveforms
%           troughAmp:   trough amplitude (uV)
%           rise2trough: rise to trough time (ms)
%           trough2peak: trough to peak time (ms)
%           FWHM:        full width at half maximum (ms)
%           isoDist:     isolation distance (well isolated, >20)
%           isiIndex:    ISI index (well isolated, <0.2)
%           Lratio;      L-ratio
%           spkNum:      number of spikes
%           meanRate:    mean firing rate (Hz)
%           isi:         histogram of isi, with filed of cnt and t (in ms)
%           acg:         ACG of spike time, with filed of cnt and t (in ms)
%       t: time vecor for plotting spike waveforms (ms)
%
% References
%   Isolation distance: Schmitzer-Torbert et al., Neuroscience, 2005
%                       Harris et al., Neuron, 2001
%   L ratio:            Schmitzer-Torbert et al., Neuroscience, 2005
%                       Schmitzer-Torbert & Redish, J Neurophysiol, 2004
%   ISI index:          Fee et al, J Neurosci Methods, 1996
%
% Takuma Kitanishi, OCU, 2017
% Modified by Hiro Miyawaki at OCU, 2018 April

function spkStats(basepath,chanMap,varargin)

% % % test path
% basepath = 'G:\tk0056\171015PhySorted'
% load('D:\data\rec\config\Kilosort\Buzsaki256.mat','chanMap0indMatrix')

%get fullpath
[~,fInfo]=fileattrib(basepath);
basepath=fInfo.Name;
params.basepath=basepath;
params.chanMap=chanMap;

% setting
fprintf('%s%s\n','Gathering spike waveforms: ',basepath)

params.dat2uV  = 5*10^6/400/(2^16/2);       % dat to micro volt conversion
% tmpfile = fullfile(basepath,'spkStatsTmp.mat');
params.nspkMax=100; %number of spikes to determine on which channel the cluster is
params.baselineMethod='mean';
params.excludeNoise=true;
params.figSaveDir=basepath;
params.matSaveDir=basepath;
params.datPath=basepath;
params.outputFigure=true; 
params.outputOnlyGood=true; 

paramList=fieldnames(params);
if mod(length(varargin),2)==1
    error('parameters must be pairs of name and value')
end
for n=1:length(varargin)/2
    name=varargin{n*2-1};
    val=varargin{n*2};
    
    ii=find(strcmpi(name,paramList));
    if isempty(ii)
        error('wrong parameter name: %s',name)
    end
    params.(paramList{ii})=val;    
end

for n=1:length(paramList)
    eval([paramList{n} '=params.' paramList{n} ';'])
end

if ~exist(figSaveDir,'dir') && outputFigure
    mkdir(figSaveDir)
end

if ~exist(matSaveDir,'dir')
    mkdir(matSaveDir)
end


% n1 = size(chanMap0indMatrix,1);     % number of channels in a shank
% n2 = size(chanMap0indMatrix,2);     % number of shanks
if iscell(chanMap)
    %cell is used for non-rectangular map (eg Buz64sp)
    chanmapLinear=[chanMap{:}]+1;
else
    chanmapLinear = chanMap(:)+1;
    temp={};
    if size(chanMap,1)==1
        chanMap=chanMap';
    end
    for n=1:size(chanMap,2)
        temp{n}=chanMap(:,n)';
    end
    chanMap=temp;
    clear temp
end
n2=length(chanMap);
% load spike sample (ss) and clusters (clu)
% [ts,clu,cluOri,Fs]=loadSpk(basepath);
% ss = round(ts*Fs);
if ~exist(fullfile(basepath,'params.py'),'file');
    error('params.py is not in %s', besepath);
end
if ~exist(fullfile(basepath,'spike_times.npy'),'file');
    error('spike_times.npy is not in %s', besepath);
end
if ~exist(fullfile(basepath,'spike_clusters.npy'),'file');
    error('spike_clusters.npy is not in %s', besepath);
end
if ~exist(fullfile(basepath,'cluster_groups.csv'),'file');
    error('cluster_groups.csv is not in %s', besepath);
end

% total recorded channel number & samplerates
paramPhy = loadParamsPy(basepath);
NchanTOT = paramPhy.n_channels_dat;
Fs=paramPhy.sample_rate;

% find a dat file
if strcmpi(params.datPath(end-3:end),'.dat') && exist(params.datPath,'file')
    datfile = dir(params.datPath);
else    
    datfile = dir(fullfile(params.datPath,'*.dat'));
    if length(datfile)~=1
        if exist(fullfile(basepath,paramPhy.dat_path),'file')
            datfile=dir(fullfile(basepath,paramPhy.dat_path));
        else
            [datFilename, datPathname] = uigetfile('*.dat', 'Select a .dat file');
            if isequal(filename,0)
                error('Dat file was not selected')
            else
                datfile = dir(fullfile(datPathname,datFilename));
            end
        end
    end
end
params.datPath=fullfile(datfile.folder,datfile.name);
bytes   = datfile.bytes;
sessName = datfile.name(1:end-4);

% spike samples
fprintf('%s loading %s\n',datestr(now),'spike_times.npy')
ss  = npy2mat(fullfile(basepath,'spike_times.npy'));
% load cluster number
fprintf('%s loading %s\n',datestr(now),'spike_clusters.npy')
cluOri = npy2mat(fullfile(basepath,'spike_clusters.npy'))';
[~,~,clu]=unique(cluOri);

% load clustering result
fprintf('%s loading %s\n',datestr(now),'cluster_groups.csv')
fid = fopen(fullfile(basepath,'cluster_groups.csv'),'r');
csv = textscan(fid,'%s %s');
fclose(fid);
id = cellfun(@str2num,csv{1,1}(2:end));
group = csv{1,2}(2:end);

% recording length (sec) calculated from the dat file size
sessLength =  bytes / (2*NchanTOT*Fs);

% setting for fread
NT1 = round((0.6/1000*Fs));  % '0.6/1000' indicates 0.6ms before spikes
NT2 = round((1.4/1000*Fs));  % '1.4/1000' indicates 1.4ms after spikes
NT  = NT1+NT2+1;             % sample number
t   = (-NT1:NT2)/Fs*1000;    % time vector (ms)

% gather spike waveforms
disp('Gathering spike waveforms')
disp(['from ' fullfile(datfile.folder,datfile.name)])
nclu = max(clu);
stats = struct;

% fid = fopen(fullfile(basepath,datfile),'r');
nSampleDat=bytes / (2*NchanTOT);
dat=memmapfile(fullfile(datfile.folder,datfile.name),'format',{'int16',[NchanTOT,nSampleDat],'raw'});


%% determine the largest amp shank by getting sparsely sampled 100 spikes
fprintf('%s Getting shank of each cluster ... ', datestr(now));
fprintf('\n  ')
prog=sprintf('%s cluster %d / %d', datestr(now),1,nclu);
fprintf(prog)
for ii = 1:nclu
    % spike timing sample of the cluster
    if mod(ii,10)==0
        fprintf(repmat('\b', 1, numel(prog)))
        prog=sprintf('%s cluster %d / %d', datestr(now),ii,nclu);
        fprintf(prog)
        if mod(ii,100)==0
            fprintf('\n  ')
            fprintf(prog)
        end
    end
    
    ssTmp = ss(clu==ii);
    n = length(ssTmp);
    
    if n<nspkMax
        nspk = 1:n;
    else
        % nspk = [1:50 n-49:n];                     % initial and last 50 spikes
        % nspk = 1:floor(n/nspkMax):floor(n/100)*99+1;    % sparse sampling
        nspk=round((0.5:nspkMax-0.5)*n/nspkMax); % avoid first and last spikes, they may have something weird
    end
    %     wavesTmp = single(nan(NT,NchanTOT,length(nspk)));
    
    %     for jj=1:length(nspk)
    % get recording data around each spike
    %         fseek(fid, 2*NchanTOT*(ssTmp(nspk(jj))-NT1), 'bof');
    %         dat = fread(fid, [NchanTOT NT], '*int16');
    %         dat = single(dat');
    wavesTmp=permute(single(reshape((dat.Data.raw(chanmapLinear,(ssTmp(nspk))+(-NT1:NT2))),[],length(nspk),NT)),[3,1,2]);
    
    if strcmpi(baselineMethod,'interpolate')
        baseline=interp1([-NT1,NT2],...
            (cat(1,mean(wavesTmp(1:3,:,:),1),mean(wavesTmp(end-2:end,:,:),1))),...
            -NT1:NT2);
    else %mean
        baseline=repmat(mean(wavesTmp([1,end],:,:),1),size(wavesTmp,1),1,1);
    end
    wavesTmp=wavesTmp-baseline;
    
    %         if size(dat,1)==NT
    %             % baseline correction
    %             wavesTmp(:,:,jj) = dat-mean([dat(1,:);dat(end,:)]);
    %         else
    %             wavesTmp(:,:,jj) = nan(NT,NchanTOT);
    %         end
    %     end
    % trim/sort channels
    %     wavesTmp = wavesTmp(:,chanmapLinear,:);
    
    % channel number with maximum amplitude
    meanWav = nanmean(wavesTmp,3);
    %     [~,maxAmpCh] = max(max(meanWav) - meanWav(NT1,:));
    [~,maxAmpCh]=max(range(meanWav)); % it's better when spikes can be positive
    
    % shank & channel(from top to bottom) number with maximum amplitude
    %     stats(ii).maxSh = ceil(maxAmpCh/n1);
    %     stats(ii).maxCh = maxAmpCh-(stats(ii).maxSh-1)*n1;
    stats(ii).maxSh = find(cellfun(@(x) ismember(chanmapLinear(maxAmpCh)-1,x),chanMap));
    stats(ii).maxCh = find(chanMap{stats(ii).maxSh}==chanmapLinear(maxAmpCh)-1);
    stats(ii).maxChRaw=chanmapLinear(maxAmpCh);
    stats(ii).group=group{ii};
end

fprintf(repmat('\b', 1, numel(prog)))
prog=sprintf('%s cluster %d / %d', datestr(now),ii,nclu);
fprintf(prog)
fprintf('\n')


%% gather spikes for each shank

% if exist(tmpfile,'file')
%     disp('Loading previously stored spkStatsTmp.mat...')
%     load(tmpfile,'shank');
% else
maxSh = [stats.maxSh];

clear shank
for ii=1:n2
    % cluster number in this shank
    if excludeNoise
        cluList = find(maxSh==ii & ~strcmpi('noise',{stats.group}));
    else
        cluList = find(maxSh==ii);
    end
    cluTmp = clu(any(clu==cluList,2));
    ssTmp  = ss(any(clu==cluList,2));
    n = length(ssTmp);
    fprintf('%s Shank: %u #spikes: %u\n',datestr(now),ii,n);
    
    if n>1
        % channel range
        chrange=chanMap{ii}+1;
        %             chrange = (ii-1)*n1+1:ii*n1;
        n1=length(chrange);
        waves=zeros(NT,n1,n);
        fprintf('   %s loading spikes...',datestr(now))
        prog='';
        fprintf('\n   ')
        for jj=1:ceil(n/1000)
            fprintf(repmat('\b', 1, numel(prog)))
            prog=sprintf('%s spike %u / %u', datestr(now),jj*1000,n);
            fprintf(prog)
            if mod(jj,100)==0
                fprintf('\n   ')
                fprintf(prog)
            end
            
            target=1+(jj-1)*1000:min(n,jj*1000);
            % get all spikes in the max amplitude shank
            %                 temp=permute(single(reshape((dat.Data.raw(chrange,(ssTmp(target))+(-NT1:NT2))),[],length(target),NT)),[3,1,2]);
            temp=permute(single(reshape((dat.Data.raw(chrange,(ssTmp(target))+(-NT1:NT2))),[],length(target),NT)),[3,1,2]);
            
            if strcmpi(baselineMethod,'interpolate')
                baseline=interp1([-NT1,NT2],...
                    (cat(1,mean(temp(1:3,:,:),1),mean(temp(end-2:end,:,:),1))),...
                    -NT1:NT2);
            else %mean
                baseline=repmat(mean(temp([1,end],:,:),1),size(temp,1),1,1);
            end
            waves(:,:,target)=temp-baseline;
            %             waves = single(nan(NT,n1,n));
        end
        fprintf(repmat('\b', 1, numel(prog)));
        prog=sprintf('%s spike %u / %u', datestr(now),target(end),n);
        fprintf(prog);
        fprintf('\n');
        
        %             for jj=1:n
        %                 % get waveform around each spike
        %                 fseek(fid, 2*NchanTOT*(ssTmp(jj)-NT1), 'bof');
        %                 dat = fread(fid, [NchanTOT NT], '*int16');
        %                 dat = single(dat');
        %
        %                 if size(dat,1) == NT
        %                     % baseline correction
        %                     dat = dat-mean([dat(1,:);dat(end,:)]);
        %                 else
        %                     dat = nan(NT,NchanTOT);
        %                 end
        %                 % trim/sort channels
        %                 dat = dat(:,chanmapLinear);
        %                 waves(:,:,jj) = dat(:,chrange);
        %             end
        
        % features for isolation distance and L ratio
        %   in Harris et al., Neuron, 2001
        %       top 3 PCs
        %   in Schmitzer-Torbert & Redish, J Neurophysiol, 2004
        %       energy and top 3 PCs on energy normalized waveform
        %   in Schmitzer-Torbert et al., Neuroscience, 2005
        %       energy and top 1 PC on energy nomalized waveform
        %
        % Here we use the method in Schmitzer-Torbert et al., Neuroscience, 2005
                       
        fprintf('%s PCA...\n',datestr(now))
        %             fet = nan(n,n1*2);
        fet = zeros(n,n1*2);
        % power of spike waveform
        e = sum(waves.^2,1);
        fet(:,1:n1) = squeeze(e)';
        % normaized waveforms
        nwaves = waves./repmat(e,NT,1,1);
        % perform pca and take the 1st component
        for jj=1:n1
            fprintf('  %s on ch %u / %u\n', datestr(now),jj,n1);
            [~,score] = pca( squeeze(permute(nwaves(:,jj,:),[2 3 1])) );
            fet(:,n1+jj) = score(:,1);
        end
        
        waves = waves * dat2uV;
        for cluIdx=1:length(cluList)
            targetClu=cluList(cluIdx);
            
            if ~outputFigure || (outputOnlyGood && ~strcmpi(stats(targetClu).group,'good'))
                doOutput=false;
            else
                doOutput=true;
            end
            %             shank(ii).waves = waves * dat2uV;   % conversion to uV
            %             shank(ii).clu = cluTmp;
            %             shank(ii).fet = fet;
            %         else
            %             shank(ii).waves = nan;
            %             shank(ii).clu   = nan;
            %             shank(ii).fet   = nan;
            %         end
            %     end
            
            %     clear waves nwaves cluTmp fet ssTmp score;
            
            %     disp('Saving tmp file...')
            %     save(tmpfile,'shank','-v7.3');
            % end
            
            % fclose(fid);
            
            
            %% mean, std, sem of spike waveform of each cluster
            % for ii=1:nclu
            %     maxSh = stats(ii).maxSh;
            %     waves = shank(maxSh).waves(:,:,shank(maxSh).clu==ii);
            %     spkNum = size(waves,3);
            spkNum = sum(cluTmp==targetClu);
            stats(targetClu).waveMean  = nanmean(waves(:,:,cluTmp==targetClu),3);
            stats(targetClu).waveStd   = nanstd(waves(:,:,cluTmp==targetClu),[],3);
            stats(targetClu).waveSem   = stats(targetClu).waveStd / sqrt(spkNum);
            %     stats(ii).waveMean  = nanmean(waves,3);
            %     stats(ii).waveStd   = nanstd(waves,[],3);
            %     stats(ii).waveSem   = stats(ii).waveStd / sqrt(spkNum);
            
            % plot waveforms for all channels in the shank
            if doOutput
                figure(1); clf
                subplot(121)

                x = stats(targetClu).waveMean;
                e = stats(targetClu).waveStd;
                pad = repmat(-1*[1:n1],NT,1);

                x2 = x(:,stats(targetClu).maxCh);
                e2 = e(:,stats(targetClu).maxCh);
                pad2 = -stats(targetClu).maxCh;

                sc = 1/100;
                hold on
                plot(t,x*sc+pad,'k',t,x2*sc+pad2,'r')
                plot(t,(x2+e2)*sc+pad2,'r:',t,(x2-e2)*sc+pad2,'r:');
                box off
                xlabel('Time (ms)')
                ylabel('Channels (channel spacing = 100uV)')
                id = cluOri(clu==targetClu);
                id = id(1);
                title(sprintf('%s%s%u%s%u%s%u%s%u%s',sessName,'_',targetClu,' (id',id,') (sh',stats(targetClu).maxSh,' ch',stats(targetClu).maxCh,')'),'Interpreter','none')
                xlim([-NT1,NT2]/Fs*1000)
                ax=fixAxis;
                plot([0 0],ax(3:4),':');

                switch lower(stats(targetClu).group)
                    case 'good'
                        textCol=[0,0.85,0];
                    case 'mua'
                        textCol=[0,0,1];
                    case 'noise'
                        textCol=[1,0,0];
                    otherwise
                        textCol=[0,0,0];
                end
                text2(-0.2,1.05,stats(targetClu).group,ax,...
                        'horizontalAlign','left','color',textCol,'fontsize',15)
            end
            % calculate waveform params
            tq = (-NT1:0.1:NT2)/Fs*1000;
            xq = interp1(t,x2,tq,'spline');  % 10-fold upsampling
            
            % flip if the spike is positive
            if abs(min(xq))<max(xq)
                xq=-xq;
            end
            
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
            %     nspk = size(waves,3);
            %     meanRate = nspk/sessLength;
            meanRate = spkNum/sessLength;
            % isolation distance
            inIdx =cluTmp==targetClu;
            %     inIdx = shank(maxSh).clu==ii;
            %     md = mahal(shank(maxSh).fet(~inIdx,:),shank(maxSh).fet(inIdx,:));
            if sum(inIdx)>size(fet,2)
                md = mahal(fet,fet(inIdx,:));
                nIn = sum(inIdx);
                nOut= sum(~inIdx);
                if nIn>nOut
                    isoDist = nan;
                else
                    sortedmd = sort(md((~inIdx)));
                    isoDist = sortedmd(sum(inIdx));
                end

                %L ratio
                df = size(fet,2);
                Lratio=sum(1-chi2cdf(md(~inIdx).^2,df))/sum(inIdx);
            else
                isoDist=nan;
                Lratio=nan;
            end
            
            % isi index (Fee et al., 1996, J Neurosci Methods)
            %     tsTmp = ts(clu==ii);
            tsTmp=ssTmp(cluTmp==targetClu)/Fs;
            isi = diff(tsTmp) *1000;  % (ms)
            isiIndex = (8/2)*( sum(isi<2)/sum(2<=isi & isi<10) );
            [isiCnt,isiBin]=hist(isi(isi<50),50);
            
            % plot the largest amplitude waveform
            if doOutput
                subplot(222)
                plot(tq,xq,'-'); hold on
                plot(tq(riseIdx),rise,'ro',tq(troughIdx),trough,'ro',tq(peakIdx),peak,'ro')
                plot([tq(FWHM1idx) tq(FWHM2idx)],[trough/2 trough/2],'r')
                hold off
                xlabel('Time (ms)')
                ylabel('Amplitude (uV)')
                title(sprintf('%s%.2f%s%.2f%s%.2f%s%.1f','Rise2tr=',rise2trough,' Tr2pk=',trough2peak,' FWHM=',FWHM,' Rate=',meanRate))

                xlim([-NT1,NT2]/Fs*1000)
                fixAxis;
                box off

                subplot(224)
                bar(isiBin,isiCnt,1)
                xlabel('ISI (ms)');
                ylabel('# spikes')
                title(sprintf('ISIidx=%.3f, IsoDist=%.1f, L_{ratio}=%.3f',isiIndex,isoDist,Lratio))
                box off
                drawnow;

                % save image
                pngImage = sprintf('%s%u%s','spkStats_',targetClu,'.png');
                epsImage = sprintf('%s%u%s','spkStats_',targetClu,'.eps');
                print('-dpng',fullfile(figSaveDir,pngImage));
                print('-depsc','-tiff',fullfile(figSaveDir,epsImage));
            end
            isiBin(end)=[];
            isiCnt(end)=[];
            [acgCnt,acgBin]=CCG(ssTmp(cluTmp==targetClu),1,1e-3*Fs, 30, Fs);
            
            % register
            stats(targetClu).id          = id;              % cluster id in Phy
            stats(targetClu).troughAmp   = -trough;         % trough amplitude (uV)
            stats(targetClu).rise2trough = rise2trough;     % rise to trough (ms)
            stats(targetClu).trough2peak = trough2peak;     % trough to peak (ms)
            stats(targetClu).FWHM        = FWHM;            % full width at half max (ms)
            stats(targetClu).isoDist     = isoDist;         % isolation distance
            stats(targetClu).isiIndex    = isiIndex;        % ISI index
            stats(targetClu).Lratio     = Lratio;           % L-ratio
            stats(targetClu).meanRate    = meanRate;        % mean firing rate (Hz)
            stats(targetClu).spkNum      = spkNum;          % number of spikes
            stats(targetClu).isi.cnt     = isiCnt;
            stats(targetClu).isi.t     = isiBin;
            stats(targetClu).acg.cnt     = acgCnt;
            stats(targetClu).acg.t     = acgBin;

            %     stats(ii).id          = id;              % cluster id in Phy
            %     stats(ii).troughAmp   = -trough;         % trough amplitude (uV)
            %     stats(ii).rise2trough = rise2trough;     % rise to trough (ms)
            %     stats(ii).trough2peak = trough2peak;     % trough to peak (ms)
            %     stats(ii).FWHM        = FWHM;            % full width at half max (ms)
            %     stats(ii).isoDist     = isoDist;         % isolation distance
            %     stats(ii).isiIndex    = isiIndex;        % ISI index
            %     stats(ii).meanRate    = meanRate;        % mean firing rate (Hz)
            %     stats(ii).spkNum      = spkNum;          % number of spikes
            %
        end
        
        
    end
end
figure(2); clf
plotParamList={'isoDist','Lratio','isiIndex','troughAmp','rise2trough','trough2peak','FWHM','meanRate'};
labelText.isoDist='Isolation distance';
labelText.Lratio='L ratio';
labelText.isiIndex='ISI index';
labelText.troughAmp='Trough amplitude (\muV)';
labelText.rise2trough='Rise to trough (ms)';
labelText.trough2peak='Trough to peak (ms)';
labelText.FWHM='FWHM (ms)';
labelText.meanRate='Mean Rate (Hz)';

for n=1:8
    subplot(4,2,n)
    [cnt,bin]=hist([stats.(plotParamList{n})],40);  
    gcnt=hist([stats(strcmpi({stats.group},'good')).(plotParamList{n})],bin);
    bar(bin,cnt,1,'k')
    hold on
    bar(bin,gcnt,1,'facecolor',[0,0.9,0])
    xlabel(labelText.(plotParamList{n}))
    box off
end
% subplot(421); hist([stats.isoDist],40);      xlabel('Isolation distance'); title(sessName)
% subplot(422); hist([stats.isiIndex],40);     xlabel('ISI index')
% subplot(423); hist([stats.troughAmp],40);    xlabel('Trough amplitude (uV)')
% subplot(424); hist([stats.rise2trough],40);  xlabel('Rise to trough (ms)')
% subplot(425); hist([stats.trough2peak],40);  xlabel('Trough to peak (ms)')
% subplot(426); hist([stats.FWHM],40);         xlabel('FWHM (ms)')
% subplot(427); hist([stats.meanRate],40);     xlabel('Mean Rate (Hz)')
print('-dpng',fullfile(figSaveDir,'spkStats_s1'));

figure(3); clf
pairList={'trough2peak','FWHM'
          'rise2trough','FWHM'
          'trough2peak','meanRate'
          'FWHM','meanRate'}
      
for n=1:4
    subplot(2,2,n)
    scatter([stats(~strcmpi({stats.group},'good')).(pairList{n,1})],...
            [stats(~strcmpi({stats.group},'good')).(pairList{n,2})],...
            36,'k')
    hold on
    scatter([stats(strcmpi({stats.group},'good')).(pairList{n,1})],...
            [stats(strcmpi({stats.group},'good')).(pairList{n,2})],...
            36,[0,0.9,0])
    xlabel(labelText.(pairList{n,1}))
    ylabel(labelText.(pairList{n,1}))
end
% subplot(221); scatter([stats.trough2peak],[stats.FWHM]); xlabel('Trough to peak (ms)'); ylabel('FWHM (ms)'); title(sessName)
% subplot(222); scatter([stats.rise2trough],[stats.FWHM]); xlabel('Rise to trough (ms)'); ylabel('FWHM (ms)'); title(sessName)
% subplot(223); scatter([stats.trough2peak],[stats.meanRate]); xlabel('Trough to peak (ms)'); ylabel('Mean rate (Hz)')
% subplot(224); scatter([stats.FWHM],[stats.meanRate]); xlabel('FWHM (ms)'); ylabel('Mean rate (Hz)')
print('-dpng',fullfile(figSaveDir,'spkStats_s2'));

% unit number map
%     unitNmap = zeros(n1,n2);
%     for ii=1:nclu
%         unitNmap(stats(ii).maxCh,stats(ii).maxSh) = unitNmap(stats(ii).maxCh,stats(ii).maxSh)+1;
%     end
bin={1:max(cellfun(@length,chanMap)),1:length(chanMap)};
allUnitNmap=hist3([[stats.maxCh];[stats.maxSh]]',bin);
[goodUnitNmap,bin]=hist3([[stats(strcmpi({stats.group},'good')).maxCh];
                          [stats(strcmpi({stats.group},'good')).maxSh]]',bin);
    

    figure(4); clf
for n=1:2
    switch n
        case 1
            unitNmap = allUnitNmap;
            utype='';
        case 2
            unitNmap = goodUnitNmap;
            utype=' good'
    end
    subplot(2,1,n)
    imagesc(unitNmap); colorbar
    pbaspect([2400 1600 1]) % probe size
%     title(sprintf('%s%s%s%u%s','Unit # map: ',sessName,'  (Total' ' unit #: ',nclu,')'));
    title(sprintf('Unit # map: %s  (Total%s unit #:%u)',sessName, utype ,sum(unitNmap(:))));
    xlabel('Shank'); ylabel('Channel')
end
    print('-dpng',fullfile(figSaveDir,'spkStats_s3'));
    print('-depsc','-tiff',fullfile(figSaveDir,'spkStats_s3'));

% save
params.generatorname=mfilename
params.updated=date;

save(fullfile(matSaveDir,'spkStats.mat'),'stats','t','sessName','chanMap','params')

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
%
% function [ts,clu,cluOri,Fs]=loadSpk(basepath,varargin)
%
% % No session specified
% if nargin==1
%     sessFlag = false;
%     disp('No session specified. Load the entire period of spike timing.')
% end
% % Session specified
% if nargin==2
%     sessFlag = true;
%     sess = varargin{1};
%
%     % load info file
%     [sessName,sessLength] = loadInfo(basepath);
%
%     % find session number
%     if isnumeric(sess)
%         sessNum = sess;
%     elseif ischar(sess)
%         sessNum = find(strcmp(sessName,sess));
%     else
%         error('Unknown input format!')
%     end
%
%     % session start/end time(sess)
%     if isempty(sessNum)
%         error('Cannot find the valid session!')
%     end
%     if sessNum==1
%         tsStart = 0;
%         tsEnd = sessLength(1);
%     else
%         cum = cumsum(sessLength);
%         tsStart = cum(sessNum-1);
%         tsEnd = cum(sessNum);
%     end
%     fprintf('%s%s%s%.1f%s%.1f%s\n','Session: ',sessName{sessNum},' (',tsStart,'-',tsEnd,' sec)')
% end
% if nargin>2
%     error('Too many inputs!')
% end
%
% % check whether previously loaded mat file exists
% matFlag = false;
% filename = fullfile(basepath,'loadSpk.mat');
% if exist(filename,'file')
%     % check timestamp of
%     load(filename,'datenum1','datenum2')
%
%     info1 = dir(fullfile(basepath,'spike_clusters.npy'));
%     info2 = dir(fullfile(basepath,'cluster_groups.csv'));
%     if datenum1 == info1.datenum && datenum2 == info2.datenum
%         matFlag = true;
%     end
% end
%
% % load spikes
% if matFlag
%     disp('Loading spikes from previously saved loadSpk.mat file...')
%     load(filename,'ts','clu','cluOri','Fs');
% else
%     disp('Loading spikes from npy/csv files...')
%
%     % get sampling rate from params.py
%     params = loadParamsPy(basepath);
%     Fs = params.sample_rate;
%
%     % load spike timing (unit in samples)
%     ssTmp  = npy2mat(fullfile(basepath,'spike_times.npy'));
%     % load cluster number
%     cluTmp = npy2mat(fullfile(basepath,'spike_clusters.npy'))';
%
%     % load clustering result
%     % 0=noise, 1=MUA, 2=Good, 3=unsorted
%     fid = fopen(fullfile(basepath,'cluster_groups.csv'),'r');
%     csv = textscan(fid,'%s %s');
%     fclose(fid);
%
%
%     id = cellfun(@str2num,csv{1,1}(2:end));
%     groupTmp = csv{1,2}(2:end);
%
%     n = length(groupTmp);
%     group = nan(n,1);
%     for ii=1:n
%         switch groupTmp{ii}
%             case 'noise'
%                 group(ii) = 0;
%             case 'mua'
%                 group(ii) = 1;
%             case 'good'
%                 group(ii) = 2;
%             case 'unsorted'
%                 group(ii) = 3;
%         end
%     end
%
%     % take only 'good' clusters
%
%     % good cluster id
%     id_good = id(group==2);
%
%     % spike index of 'good' clusters
%     idx = logical(sum(cluTmp(:) == id_good',2));
%
%     ts  = ssTmp(idx)/Fs;
%     cluOri = cluTmp(idx);
%
%     % sort cluster number
%     [~,~,clu]=unique(cluOri);
%
%     % date modified
%     info1 = dir(fullfile(basepath,'spike_clusters.npy'));
%     datenum1 = info1.datenum;
%     info2 = dir(fullfile(basepath,'cluster_groups.csv'));
%     datenum2 = info2.datenum;
%
%     % save
%     save(fullfile(basepath,'loadSpk.mat'),'ts','clu','cluOri','Fs','datenum1','datenum2')
% end
%
% % select ts in the session
% if sessFlag
%     idx = tsStart<ts & ts<tsEnd;
%     ts = ts(idx)-tsStart;
%     clu = clu(idx);
%     cluOri = cluOri(idx);
% end
%
%
% % read other parameters
% % amp            = npy2mat(fullfile(path,'amplitudes.npy'));
% % chanMap        = npy2mat(fullfile(path,'channel_map.npy'));
% % chanPos        = npy2mat(fullfile(path,'channel_positions.npy'));
% % pcFetInd       = npy2mat(fullfile(path,'pc_feature_ind.npy'));
% % pcFet          = npy2mat(fullfile(path,'pc_features.npy'));
% % simTemplates   = npy2mat(fullfile(path,'similar_templates.npy'));
% % spkTemplates   = npy2mat(fullfile(path,'spike_templates.npy'));
% % templateFetInd = npy2mat(fullfile(path,'template_feature_ind.npy'));
% % templateFet    = npy2mat(fullfile(path,'template_features.npy'));
% % template       = npy2mat(fullfile(path,'templates.npy'));
% % templateInd    = npy2mat(fullfile(path,'templates_ind.npy'));
% % templateUnw    = npy2mat(fullfile(path,'templates_unw.npy'));
% % whitening      = npy2mat(fullfile(path,'whitening_mat.npy'));

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

%%
% function h=text2(Xpos,Ypos,String,Axis,varargin)
% wrapper of text to specify text as in relative position
%
% Hiro Miyawaki at UWM, 2011-
%
function h=text2(Xpos,Ypos,String,Axis,varargin)

    if nargin<4
        Axis=axis();
    end
    
    % just for backward compatibility
    if length(varargin)==1 && iscell(varargin)
        Property=varargin{1};
    else
        Property=varargin;
    end

    H=text(Axis(1:2)*[1-Xpos;Xpos],Axis(3:4)*[1-Ypos;Ypos],String,Property{:});
    if nargout>0
        h=H;
    end
%%    
function ax=fixAxis
    ax=axis();
    xlim(ax(1:2));
    ylim(ax(3:4));
    
%%
function t = date
%DATE   Current date as character vector.
%   S = DATE returns a character vector containing the date in dd-mmm-yyyy format.
%
%   See also NOW, CLOCK, DATENUM.

%   Copyright 1984-2016 The MathWorks, Inc.

c = clock;
mths = ['Jan';'Feb';'Mar';'Apr';'May';'Jun';'Jul';
        'Aug';'Sep';'Oct';'Nov';'Dec'];
d = sprintf('%.0f',c(3)+100);
t = [d(2:3) '-' mths(c(2),:) '-' sprintf('%.0f',c(1))];
    