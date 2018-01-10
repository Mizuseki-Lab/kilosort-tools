clear
clc

datDir='C:/Users/Teppan/Documents/Miyawaki/data/magician/removeNoise/';

datFile=dir(fullfile(datDir,'*.dat'));
nCh=112;
nSample=datFile.bytes/2/nCh;

lfp=memmapfile(fullfile(datFile.folder,datFile.name), 'Format',{'int16',[nCh,nSample],'raw'});

nFil=33;
nBef=16;
nAft=24;
for n=1:4
    chMap{n}=(1:8)+8*(n-1);
end
for n=5:10
    chMap{n}=(1:10)+32+10*(n-5);
end

usedCh=[chMap{:}];
nUsedCh=length(usedCh);
%%

nFet=3;

%%
kwikDir='C:/Users/Teppan/Documents/Miyawaki/data/magician/spikeDetekt/';
kiloDir='C:/Users/Teppan/Documents/Miyawaki/data/magician/removeNoise/';

kwikfile=dir(fullfile(kwikDir,'*.kwik'));
kwikfile=fullfile(kwikfile.folder,kwikfile.name);
%%
for type=1:2
    
    if type==1
        fprintf('%s loading npy files\n',datestr(now))

        allRes=npy2mat(fullfile(kiloDir,'spike_times.npy'));
        allClu=npy2mat(fullfile(kiloDir,'spike_clusters.npy'));

        allCluList=readtable(fullfile(kiloDir,'cluster_groups.csv'));
        allCluList=allCluList.cluster_id(cellfun(@(x) ~strcmpi(x,'noise'),allCluList.group));
        nClu=length(allCluList);
        fprintf('%s detecting max amp ch for each cluster\n',datestr(now))
        parfor n=1:nClu
            cID=allCluList(n);
            sIdx=find(allClu==cID);
            sIdx(sIdx<nBef+nFil | sIdx>nSample-nAft-nFil)=[];
            sIdx=sort(sIdx(randperm(length(sIdx),min(length(sIdx),3000))));

            wave=single(reshape(lfp.Data.raw(usedCh,allRes(sIdx)+(-nBef-nFil:nAft+nFil)),nUsedCh,length(sIdx),nBef+nAft+1+2*nFil));
            wave=wave-medfilt1(single(wave),2*nFil+1,[],3);
            wave=(squeeze(mean(wave,2))');
            [~,maxCh(n)]=max(max(wave)-min(wave));
        end
    end
    
    for sh=1:length(chMap)
        fprintf('%s proccesing on shank %d\n',datestr(now),sh)
        
        if type==1
            rootDir=kiloDir;
            nameCore=fullfile(rootDir,'magician-kilo');
            
            cluList=allCluList(ismember(maxCh,chMap{sh}));
            res=allRes(ismember(allClu,cluList));
            clu=allClu(ismember(allClu,cluList));
        elseif type==2
            rootDir=kwikDir;
            nameCore=fullfile(rootDir,'magician-kwik');
            
            clu=h5read(kwikfile,sprintf('/channel_groups/%d/spikes/clusters/main',sh-1));
            res=h5read(kwikfile,sprintf('/channel_groups/%d/spikes/time_samples',sh-1));
            
            cluList=unique(clu);
            cluGrp={};
            for cIdx=1:length(cluList)
                temp=h5readatt(kwikfile,sprintf('/channel_groups/%d/clusters/main/%d/',sh-1,cluList(cIdx)),'cluster_group');
                cluGrp(cIdx)=h5readatt(kwikfile,sprintf('/channel_groups/%d/cluster_groups/main/%d/',sh-1,temp),'name');
            end
            cluList=cluList(~strcmpi(cluGrp,'noise'));
            res=double(res(ismember(clu,cluList)));
            clu=clu(ismember(clu,cluList));
            
        else
            continue
        end
        
        clu(res<nFil+nBef | res>nSample-nFil-nAft)=[];
        res(res<nFil+nBef | res>nSample-nFil-nAft)=[];
        
        nSpk=length(res);
        nCh=length(chMap{sh});
        fprintf('    %s loading and filtering spikes (%d spikes)\n',datestr(now),nSpk)
        spk=zeros(nCh,nSpk,nBef+nAft+1);
        parfor n=1:nSpk
            temp=lfp.Data.raw(chMap{sh},res(n)+(-nBef-nFil:nAft+nFil));
            temp=temp-int16(medfilt1(single(temp),nFil*2+1,[],2));
            spk(:,n,:)=temp(:,nFil+1:end-nFil)
        end
        
        fprintf('    %s getting pca score\n',datestr(now));
        fet=zeros(nSpk,nFet*nCh);
        for ch=1:nCh
            temp=single(squeeze(spk(ch,:,:)));
            [u,s,v]=svd(temp-mean(temp),'econ');
            s=diag(s);
            fet(:,(1:nFet)+nFet*(ch-1))=u(:,1:nFet).*s(1:nFet)';
        end
        
        fprintf('    %s saving results\n',datestr(now));
        save(sprintf('%s_shank%02d.mat',nameCore,sh),'fet','clu','res','-v7.3')
        %save(sprintf('%s_shank%02d-spk.mat',nameCore,sh),'spk','-v7.3')
        
    end
end


