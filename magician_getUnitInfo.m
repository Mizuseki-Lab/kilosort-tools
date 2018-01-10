clear
clc
%%
FileInfo=dir(fullfile('/Users/miyawaki/Desktop/removeNoise/','*dat'));
NChInDat=112;
dat = memmapfile(fullfile(FileInfo.folder,FileInfo.name), 'Format', {'int16', [NChInDat, (FileInfo.bytes/NChInDat/2)], 'data'});


tBin=30:60:size(dat.Data.data,2)/20e3;
nSpk=10000;

for n=1:4
    chMap{n}=(1:8)+8*(n-1);
end
for n=5:10
    chMap{n}=(1:10)+32+10*(n-5);
end

nBef=16;
nAft=24;
nFil=33;
%%
path='/Users/miyawaki/Desktop/removeNoise/';

spk=double(npy2mat(fullfile(path,'spike_times.npy')));
clu=npy2mat(fullfile(path,'spike_clusters.npy'));

fh = fopen(fullfile(path,'cluster_groups.csv'),'r');
cluGrp = textscan(fh,'%s %s');
fclose(fh);

temp=cluGrp{1}(~strcmpi(cluGrp{2},'noise'));
cluList=cellfun(@str2num,temp(2:end));
cluGrp=cluGrp{2}(~strcmpi(cluGrp{2},'noise'))
cluGrp(1)=[];

usedCh=cat(2,chMap{:});

cnt=0;
for cIdx=1:length(cluList)
    disp([datestr(now) sprintf(' %d/%d of kilosort cluster', cIdx,length(cluList))])
    cluID=cluList(cIdx);    
    subset=spk(clu==cluID);
    
    tempAcg=CCG(subset,1,20,30,20e3);

    if ~strcmpi(cluGrp{cIdx},'good') && ...
            (mean(tempAcg(31+(-2:2)))>=mean(tempAcg) || length(subset)<2000 || ...
             mean(tempAcg)<2 || max(tempAcg)<5)
        continue
    end
    
    cnt=cnt+1;
        
    unitInfo.kilo(cnt).id=cluID;
    unitInfo.kilo(cnt).acg=tempAcg;    
    unitInfo.kilo(cnt).n=length(subset);
    unitInfo.kilo(cnt).rate=hist(subset/20e3,tBin)/60;
    unitInfo.kilo(cnt).type=lower(cluGrp{cIdx});
    
    if length(subset)>nSpk
        example=sort(subset(randperm(length(subset),nSpk)));
    else
        example=subset;
    end
    
    wave=double(dat.Data.data(usedCh,example+(-nBef-nFil:nAft+nFil)));
    wave=double(reshape(wave,size(wave,1),size(example,1),[]));
    wave=wave-medfilt1(wave,nFil*2+1,[],3);
    wave(:,:,[1:nFil,end-nFil+1:end])=[];        
      
    tempMean=squeeze(mean(wave,2));
    [unitInfo.kilo(cnt).amp,unitInfo.kilo(cnt).ch]=max(range(tempMean,2));
    
    wave=wave(chMap{cellfun(@(x) ismember(unitInfo.kilo(cnt).ch,x),chMap)},:,:);
    unitInfo.kilo(cnt).wave.mean=squeeze(mean(wave,2))';
    unitInfo.kilo(cnt).wave.std=squeeze(std(wave,0,2))';
end
%%
filename='/Users/miyawaki/Desktop/kwik/magician_2017-09-19_06-10-17.kwik';
cnt=0;

for shankIdx=1:10
    disp([datestr(now) sprintf(' %d/%d of kwik shank', shankIdx,10)])
    spk=double(h5read(filename, sprintf('/channel_groups/%d/spikes/time_samples',shankIdx-1)));
    clu=h5read(filename, sprintf('/channel_groups/%d/spikes/clusters/main',shankIdx-1));

    cluList=unique(clu);
    cluGrp={};
    for cIdx=1:length(cluList)
        temp=h5readatt(filename, sprintf('/channel_groups/%d/clusters/main/%d/',shankIdx-1,cluList(cIdx)),'cluster_group');

        cluGrp(cIdx)=h5readatt(filename,sprintf('/channel_groups/%d/cluster_groups/main/%d/',shankIdx-1,temp),'name');    
    end

    cluList=cluList(~strcmpi(cluGrp,'noise'));
    cluGrp=cluGrp(~strcmpi(cluGrp,'noise'));


    for cIdx=1:length(cluList)
        disp(['    ' datestr(now) sprintf(' %d/%d of %d shank', cIdx,length(cluList),shankIdx)])

        cluID=cluList(cIdx);
        subset=spk(clu==cluID);
        tempAcg=CCG(subset,1,20,30,20e3);

        if ~strcmpi(cluGrp{cIdx},'good') && ...
                (mean(tempAcg(31+(-2:2)))>=mean(tempAcg) || length(subset)<2000 || ...
                 mean(tempAcg)<2 || max(tempAcg)<5)
            continue
        end
        
        cnt=cnt+1;
        cluID=cluList(cIdx);
        subset=spk(clu==cluID);
        unitInfo.kwik(cnt).id=[shankIdx,cluID];
        unitInfo.kwik(cnt).acg=tempAcg;    
        unitInfo.kwik(cnt).n=length(subset);
        unitInfo.kwik(cnt).rate=hist(subset/20e3,tBin)/60;
        unitInfo.kwik(cnt).type=lower(cluGrp{cIdx});

        if length(subset)>nSpk
            example=sort(subset(randperm(length(subset),nSpk)));
        else
            example=subset;
        end

        example(example<nBef+nFil | example>size(dat.Data.data,2)-nAft-nFil)=[];
        
        wave=double(dat.Data.data(chMap{shankIdx},example+(-nBef-nFil:nAft+nFil)));
        wave=double(reshape(wave,size(wave,1),size(example,1),[]));
        wave=wave-medfilt1(wave,nFil*2+1,[],3);
        wave(:,:,[1:nFil,end-nFil+1:end])=[];

        unitInfo.kwik(cnt).wave.mean=squeeze(mean(wave,2))';
        unitInfo.kwik(cnt).wave.std=squeeze(std(wave,0,2))';

        [unitInfo.kwik(cnt).amp,idx]=max(range(unitInfo.kwik(cnt).wave.mean,1));
        unitInfo.kwik(cnt).ch=chMap{shankIdx}(idx);
    end
end
%%

save('~/data/OCU/implanted/magician/unitInfo.mat','unitInfo','-v7.3')

clear

for n=1:4
    chMap{n}=(1:8)+8*(n-1);
end
for n=5:10
    chMap{n}=(1:10)+32+10*(n-5);
end


rootDir='~/data/OCU/implanted/magician/';
load(fullfile(rootDir,'unitInfo.mat'));
for typeIdx=2
    if typeIdx==1
        type='kilo';

        temp=zeros(size(unitInfo.(type)));
        for sh=1:10
            temp=temp+sh*ismember([unitInfo.(type).ch],chMap{sh});
        end
        target=[temp',[unitInfo.(type).id]'];        
    
    elseif typeIdx==2
        type='kwik';
        
        target=cat(1,unitInfo.(type).id);    
    else
        continue
    end
    for sh=1:10
        fprintf('%s processing shank %d of %s\n', datestr(now), sh, type)
        load(fullfile(rootDir,sprintf('magician-%s_shank%02d.mat',type,sh)))
        
        subTarget=target(target(:,1)==sh,2);

        for cIdx=1:length(subTarget)
            fprintf('    %s processing cluster %d of %d\n', datestr(now), cIdx, length(subTarget))
            cID=subTarget(cIdx);

            inFet=fet(clu==cID,:);
            outFet=fet(clu~=cID,:);

            md=mahal(outFet,inFet);
            Lratio=sum(chi2cdf(md,size(inFet,2),'upper'))/size(inFet,1);

            md=sort(md,'ascend');
            if length(md)<size(inFet,1)
                isoDist=max(md);
            else
                isoDist=md(size(inFet,1));
            end
            
            uIdx=find(target(:,1)==sh & target(:,2)==cID);
            
            unitInfo.(type)(uIdx).Lratio=Lratio;
            unitInfo.(type)(uIdx).isoDist=isoDist;    
        end    
    end
end
save(fullfile(rootDir,'unitInfo.mat'),'unitInfo','-v7.3')