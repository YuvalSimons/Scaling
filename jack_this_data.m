function [data,bootidx]=jack_this_data(data,idx,inout)
if ~exist('inout','var')
    inout='in';
end

% Option to randomly choose which blocks are included and excluded instead
% of using their 10 modulus
ri=0;
if idx>=10 ri=randi(RandStream('mt19937ar','seed',idx),10,size(data.SNPs.blk)); end

if strcmp(inout,'in')
   idx=find(mod(data.SNPs.blk+ri,10)~=mod(idx,10));
else
   idx=find(mod(data.SNPs.blk+ri,10)==mod(idx,10));
end

data.SNPs=data.SNPs(idx,:);
data.traits.n=arrayfun(@(i) sum(data.SNPs.id==i),1:data.nt)';

bootidx=idx;

end

function blk=get_blocks(chr,bp)
filename="blocks.tsv";
%if ~exist(filename,'file')
%    filename='blocks.tsv';
%end
    
tmp=readmatrix(filename,'filetype','text','TrimNonNumeric',true);

blocks=array2table(tmp,'variablenames',{'chr','down','up'});

blk=zeros(size(chr));
blknum=length(blocks.chr);

for i=1:blknum
    idx=find((blocks.chr(i)==chr).*(blocks.down(i)<=bp).*(blocks.up(i)>bp));
    blk(idx)=i;%blk=[blk;idx];
end
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

%if strcmp(inout,'in')
%newid=cumsum(data.traits.n~=0).*(data.traits.n~=0);
%data.SNPs.id=newid(data.SNPs.id);
%data.nt=data.nt-sum(data.traits.n==0);
%data.traits(data.traits.n==0)=[];
%end