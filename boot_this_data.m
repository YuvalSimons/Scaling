function [data,idx]=boot_this_data(data,idx)

if ~exist('idx','var')
    seed=get_seed();
    disp(['seed ',num2str(seed,15)]);
    rng(seed);
    blk_drwn=mnrnd(1703,ones(1,1703)/1703);
    idx=[];
    for b=find(blk_drwn)
        idx=[idx;repmat(find(data.SNPs.blk==b),blk_drwn(b),1)];
    end
end

data.SNPs=data.SNPs(idx,:);
data.traits.n=arrayfun(@(i) sum(data.SNPs.id==i),1:data.nt)';

end
    

function seed=get_seed()
    [~,host]=system('hostname');
    hostnum=uint8(host);
    seed=round(mod(datenum(datetime('now'))*1e5,1)*1e5);
    for i=1:length(hostnum)
       seed=seed+round(mod(datenum(datetime('now'))*1e5,1)*1e5);
       seed=powermod(seed,hostnum(i)+round(mod(datenum(datetime('now'))*1e5,1)*10),2^32-round(mod(datenum(datetime('now'))*1e5,1)*1e3));
    end
end

