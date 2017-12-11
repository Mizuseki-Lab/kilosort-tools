function runKiloSort(basepath,chanMap,varargin)
% basepath : directory containing dat file
% chanMap : path of channel map file
% varargin : pair of option names and values such as
%   NchanTOT : total number of channels
%   Nchan : number of active channels

% check for overwriting

if exist(fullfile(basepath,'rez.mat'),'file')
    str = input('Kilosort result (rez.mat) already exists! Overwrite? (y/n)','s');
    if strcmp(str,'n') || strcmp(str,'n')
        disp('Processing stopped')
        return;
    end
end

if mod(length(varargin),2)~=0
    error('Options should be pairs of name and value'); 
end

addpath(genpath('C:\Users\Teppan\Documents\KiloSort')) % path to kilosort folder
addpath(genpath('C:\Users\Teppan\Documents\KiloSort\npy-matlab')) % path to npy-matlab scripts

ops = config(basepath,chanMap,varargin{:});

tic; % start timer
%
if ops.GPU
    gpuDevice(1); % initialize GPU (will erase any existing GPU arrays)
end


[rez, DATA, uproj] = preprocessData(ops); % preprocess data and extract spikes for initialization
rez                = fitTemplates(rez, DATA, uproj);  % fit templates iteratively
rez                = fullMPMU(rez, DATA);% extract final spike times (overlapping extraction)

save(fullfile(ops.root,  'rez.mat'), 'rez', '-v7.3');

% Remove huge amplitude noises (Takuma)
% rez = removeNoise(rez);

% save python results file for Phy
rezToPhy(rez, ops.root);

% remove temporary file
delete(ops.fproc);

function ops = config(basepath,chanMap,varargin)

% define the channel map as a filename (string) or simply an array
ops.chanMap             = chanMap; % make this file using createChannelMapFile.m
ops.criterionNoiseChannels = 0.01; % fraction of "noise" templates allowed to span all channel groups (see createChannelMapFile for more info).
% ops.chanMap = 1:ops.Nchan; % treated as linear probe if a chanMap file

ops.GPU                 = 1; % whether to run this code on an Nvidia GPU (much faster, mexGPUall first)
ops.parfor              = 0; % whether to use parfor to accelerate some parts of the algorithm
ops.verbose             = 1; % whether to print command line progress
ops.showfigures         = 1; % whether to plot figures during optimization

ops.datatype            = 'dat';  % binary ('dat', 'bin') or 'openEphys'
datfile = dir(fullfile(basepath,'*.dat'));
if length(datfile)<1; error('No dat file in the folder!'); end
if length(datfile)>1; error('Two or more dat files in the folder!'); end

ops.fbinary = fullfile(datfile.folder,datfile.name);
ops.fproc   = fullfile(datfile.folder,'temp_wh.dat');
ops.root    = fullfile(datfile.folder,'\');
disp(ops.fbinary)
% ops.fbinary             = 'C:\data\oko2\tk0016\Kilosort\170223-1\tk0016-170223-1.dat'; % will be created for 'openEphys'
% ops.fproc               = 'C:\data\oko2\tk0016\Kilosort\170223-1\temp_wh.dat'; % residual from RAM of preprocessed data
% ops.root                = 'C:\data\oko2\tk0016\Kilosort\170223-1\'; % 'openEphys' only: where raw files are

cMap=load(chanMap);
if ~isfield(cMap,'fs')
    ops.fs                  = 20000;  % sampling rate (omit if already in chanMap file)
end

ops.NchanTOT            = max(cMap.chanMap(:));     % total number of channels (omit if already in chanMap file)
ops.Nchan               = length(cMap.chanMap(:));     % number of active channels (omit if already in chanMap file)
ops.Nfilt               = ceil(ops.Nchan*6/32)*32;    % number of clusters to use (2-4 times more than Nchan, should be a multiple of 32)
ops.nNeighPC            = 6;      % visualization only (Phy): number of channnels to mask the PCs, leave empty to skip (12)
ops.nNeigh              = 8;      % visualization only (Phy): number of neighboring templates to retain projections of (16)

% options for channel whitening
ops.whitening           = 'full'; % type of whitening (default 'full', for 'noSpikes' set options for spike detection below)
ops.nSkipCov            = 1;      % compute whitening matrix from every N-th batch (1)
ops.whiteningRange      = 32;     % how many channels to whiten together (Inf for whole probe whitening, should be fine if Nchan<=32)

% other options for controlling the model and optimization
ops.Nrank               = 3;     % matrix rank of spike template model (3)
ops.nfullpasses         = 6;     % number of complete passes through data during optimization (6)
ops.maxFR               = 20000; % maximum number of spikes to extract per batch (20000)
ops.fshigh              = 300;   % frequency for high pass filtering
% ops.fslow             = 2000;  % frequency for low pass filtering (optional)
ops.ntbuff              = 64;    % samples of symmetrical buffer for whitening and spike detection
ops.scaleproc           = 200;   % int16 scaling of whitened data
ops.NT                  = 32*1024+ ops.ntbuff; % this is the batch size (try decreasing if out of memory)
% for GPU should be multiple of 32 + ntbuff

% the following options can improve/deteriorate results.
% when multiple values are provided for an option, the first two are beginning and ending anneal values,
% the third is the value used in the final pass.
ops.Th               = [6 12 12];    % threshold for detecting spikes on template-filtered data ([6 12 12])
ops.lam              = [10 30 30];    % large means amplitudes are forced around the mean ([10 30 30])
ops.nannealpasses    = 4;            % should be less than nfullpasses (4)
ops.momentum         = 1./[20 400];  % start with high momentum and anneal (1./[20 1000])
ops.shuffle_clusters = 1;            % allow merges and splits during optimization (1)
ops.mergeT           = .1;           % upper threshold for merging (.1)
ops.splitT           = .1;           % lower threshold for splitting (.1)

% options for initializing spikes from data
ops.initialize      = 'no';    %'fromData' or 'no'
ops.spkTh           = -6;      % spike threshold in standard deviations (4)
ops.loc_range       = [3  1];  % ranges to detect peaks; plus/minus in time and channel ([3 1])
ops.long_range      = [30  6]; % ranges to detect isolated peaks ([30 6])
ops.maskMaxChannels = 5;       % how many channels to mask up/down ([5])
ops.crit            = .65;     % upper criterion for discarding spike repeates (0.65)
ops.nFiltMax        = 10000;   % maximum "unique" spikes to consider (10000)

% load predefined principal components (visualization only (Phy): used for features)
dd                  = load('PCspikes2.mat'); % you might want to recompute this from your own data
ops.wPCA            = dd.Wi(:,1:7);   % PCs

% options for posthoc merges (under construction)
ops.fracse  = 0.1; % binning step along discriminant axis for posthoc merges (in units of sd)
ops.epu     = Inf;

ops.ForceMaxRAMforDat   = 120e9; % maximum RAM the algorithm will try to use; on Windows it will autodetect.

for n=1:length(varargin)/2
    if ~isfield(ops,varargin{2*n-1})
        error(['Invalid option name :' ,varargin{2*n-1}]);
    else
        ops.(varargin{2*n-1})=varargin{2*n};
    end
end


