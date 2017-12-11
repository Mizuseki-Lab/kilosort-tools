function removeNoiseCluFromRez(path,varargin)
% remove noise clusters from rez files PERMANENTLY. Use with care!
%
% removeNoiseCluFromRez(path,[options])   
% path : directory containing rez.mat and other kilosort files
% options : pairs of option names and values
%   backup : (boolean) back up files before removing
%   update : (boolean) update npy files as well

    if mod(length(varargin),2)~=0
        error('Stopped due to wrong option')
    end
    
    backupNpy=true;
    updateNpy=true;
    for idx=1:length(varargin)/2
        name=lower(varargin{idx*2-1});
        val=varargin{idx*2};
        
        switch name
            case 'backup'
                backupNpy=val;
            case 'update'
                updateNpy=val;
            otherwise
            error(['Wrong option: ' name])
        end
    end


    if backupNpy
        backupDir=fullfile(path,'backup');
        n=0;
        while exist(backupDir,'dir')
            n=n+1;
            backupDir=fullfile(path,['backup' num2str(n)]);
        end
        mkdir(backupDir);
        
        backupFiles=dir(fullfile(path,'*.npy'));
        backupFiles={'rez.mat','cluster_groups.csv',backupFiles.name};

        for fIdx=1:length(backupFiles)
            disp([datestr(now) ' backing up ' backupFiles{fIdx}])
            status=copyfile(fullfile(path,backupFiles{fIdx}),...
                            fullfile(backupDir,backupFiles{fIdx}));
             if ~status
                 error(['Can not back up ' ,backupFiles{fIdx}])
             end             
        end
    else
            disp(['Npy files were NOT backed up' backupFiles{fIdx}])
    end


    disp([datestr(now) ' Loading rez.mat'])
    load(fullfile(path,'rez.mat'));

%     disp([datestr(now) 'Loading spike times'])
%     spk=npy2mat(fullfile(path,'spike_times.npy'));
%     spk=uint64(rez.st3(:,1));

    
    disp([datestr(now) ' Loading cluster ID'])
    clu=npy2mat(fullfile(path,'spike_clusters.npy'));
    
    fh = fopen(fullfile(path,'cluster_groups.csv'),'r');
    cluGrp = textscan(fh,'%s %s');
    fclose(fh);

    noiseCluList = cellfun(@str2num,cluGrp{1}(strcmp(cluGrp{2},'noise')));

    removeSpk=ismember(clu,noiseCluList);


    rez.st3(removeSpk,:)=[];
    rez.cProj(removeSpk,:)=[];
    rez.cProjPC(removeSpk,:,:)=[];
    
    disp([datestr(now) ' Saving updated rez.mat'])
    save(fullfile(path,'rez.mat'),'rez','-v7.3');
        
    if updateNpy
        rez.st3(:,5)=clu(~removeSpk);
        disp([datestr(now) ' re-generating npy files'])
        rezToPhy(rez, path);
        
        disp([datestr(now) ' Updating cluster_groups.csv'])
        fh = fopen(fullfile(path,'cluster_groups.csv'),'w');
        for cIdx=1:length(cluGrp{1})
            if ~strcmp(cluGrp{2}(cIdx),'noise')
                fprintf(fh,'%s\t%s\n',cluGrp{1}{cIdx},cluGrp{2}{cIdx});
            end
        end
        fclose(fh);
    else
            disp(['Npy files were NOT updated. To update, use rezToPhy() provided by kilosort'])
    end
    
    disp('')
    disp([datestr(now) ' DONE!'])
