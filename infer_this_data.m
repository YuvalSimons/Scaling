function infer_this_data(datafile,fileout,model,options)

arguments
datafile char
fileout char
model char = "SSD" 
options.boot double {mustBeInteger} = 0  %bootstrap number (0 indicates no bootstrap)
options.jack double {mustBeInteger} = -1 %jackknife number (-1 indicates no jackknife)
options.eps double = 0.005  %penalty strength
options.k double {mustBeInteger} = 4    %no. of knots minus one
options.minx double = 0.01  %minimal included MAF
options.maxld double = 500  %maximal included LD score
options.minsnps double {mustBeInteger} = 100 % minimal number of SNPs per included trait
options.ex (1,:) = [] % Excluded traits
options.splinefile char = 'splines.mat' %file with precalculated normalizing splines
options.tablefile char = 'tables.mat'  %file with precalculated prob of inclusion at the LD threshold for different MAF bins
options.funevals double = 24000 %number of function evals in likelihood maximization
options.maxgrid = 0 %location of knot with largest log10(s)
options.mingrid = -8 %location of knot with smallest log10(s)
options.maxlogs = 0 %maximal included log10(s)
options.minlogs = -8 %minimal included log10(s)
options.targetsize = 3e8 %estimated size of functional genome for the penalty
options.pout = 0.001 %Proportion of outliers SNPs (set to 0 for inference without accounting for outliers)
options.zout = 100 % mean z-score of outliers
end

%grid of C values
cvec=3:0.01:9;

% grid of selection coefficients
logs=(options.minlogs:1/16:options.maxlogs)'; %(-8:1/16:-1)'; 
s=10.^logs;
ls2i=@(ls) min(81,max(1,16*(ls+6)+1));

load(datafile,'data');
data.eps=options.eps; data.k=options.k;
data.logs=logs;
data.funevals=options.funevals;
data.mingrid=options.mingrid; data.maxgrid=options.maxgrid; 
data.targetsize=options.targetsize;
data.pout=options.pout;
data.zout=options.zout;

%Temporary code to see if the inference converges
seed=get_seed();
disp(['seed ',num2str(seed,15)]);
rng(seed);

[filepath,name,ext]=fileparts(datafile);
datafile=[name,ext];



data=preprocess(data,options.minx,options.maxld,options.minsnps,options.ex);


data.shtot=readmatrix('shtot.csv');
%Expand range
data.shtot=data.shtot(ls2i(logs));
%data.shtot=[s(logs<-6)*1e6*data.shtot(1);data.shtot];
%data.shtot=[data.shtot;ones(sum(logs>-1),1)*data.shtot(end)];

% if strcmp(options.pfile,'')
%     if ~isempty(filepath) options.pfile=[filepath,'/',name,'_p_ics.mat']; else options.pfile=[name,'_p_ics.mat']; phitfile=[name,'_phit_ics.mat']; end
% end
% if strcmp(options.phitfile,'')
%     if ~isempty(filepath) options.phitfile=[filepath,'/',name,'_phit_ics.mat']; else options.phitfile=[name,'_phit_ics.mat']; end
% end
% if exist(options.pfile,'file') 
%     load(options.pfile); 
% else
%     p_ics=save_ptable(datafile,'minx',options.minx,'maxld',options.maxld,'minsnps',options.minsnps,'ex',options.ex,'splinefile',options.splinefile,'tablefile',options.tablefile);
% end
% if exist(options.phitfile,'file') 
%     load(options.phitfile); 
% else
%     phit_ics=save_phittable(datafile,'minx',options.minx,'maxld',options.maxld,'minsnps',options.minsnps,'ex',options.ex,'splinefile',options.splinefile,'tablefile',options.tablefile);
% end

str='';

if data.eps>0
    str=[str,'_eps',num2str(data.eps)];
end
if data.k~=4
    str=[str,'_k',num2str(data.k)];
end
if data.mingrid~=-6
    str=[str,'_mingrid',num2str(data.mingrid)];
end
if data.maxgrid~=-2
    str=[str,'_maxgrid',num2str(data.maxgrid)];
end
if options.minlogs~=-8
    str=[str,'_minlogs',num2str(options.minlogs)];
end
if options.maxlogs~=0
    str=[str,'_maxlogs',num2str(options.maxlogs)];
end
if options.targetsize~=3e9
    str=[str,'_targetsize',num2str(options.targetsize)];
end
if options.pout~=0
    str=[str,'_pout',num2str(options.pout)];
    str=[str,'_zout',num2str(options.zout)];
end


[folderout,filename,ext]=fileparts(fileout);
if ~strcmp(folderout,'') folderout=[folderout,'/']; mkdir(folderout); end

bootidx=[]; outidx=[]; pres=[];
if options.boot>0
    [data,bootidx]=boot_this_data(data);
    folderout=[folderout,'/boot/'];
    str=[str,'_boot',num2str(options.boot)];
    disp(['size of dataset is ',num2str(length(data.SNPs.x))]);
elseif options.jack>=0
    [dataout,outidx]=jack_this_data(data,options.jack,'out');
    [data,bootidx]=jack_this_data(data,options.jack);
    folderout=[folderout,'/jack/'];
    str=[str,'_jack',num2str(options.jack)];
    disp(['size of dataset is ',num2str(length(data.SNPs.x))]);
end

mkdir(folderout);
filename=[folderout,filename,'_',model,str,'.mat'];
disp(filename);

%We apply the bootstrap/jackknife
% if ~isempty(bootidx) 
%     p_ics=p_ics(booidx,:,:);
%     phit_ics=phit_ics(booidx,:,:);
% end
% %We stack the p tables into 2*2 matrices
% f=@(m) reshape(m,[length(data.SNPs.x)*length(cvec),length(logs)]);
% p_ics=f(p_ics);
% phit_ics=f(phit_ics);








if ~exist(filename,"file")
tic;
inf=[];
if strcmp(model,'SSD')
    %[p_ics,phit_ics]=make_ptables(data);
    [p_ics,pout_i]=make_ptable(data,'splinefile',options.splinefile,'tablefile',options.tablefile);
    phit_ics=make_phittable(data,'splinefile',options.splinefile,'tablefile',options.tablefile);

    %We need to reshape the 3d tables into 2d matrices
    f=@(m) reshape(m,[length(data.SNPs.x)*length(cvec),length(logs)]);
    p_ics=f(p_ics);
    phit_ics=f(phit_ics);

    
    inf=do_inference(data,p_ics,phit_ics,pout_i);    
    

elseif strcmp(model,'TSD')
    for trt=1:data.nt
        tmp=subset_traits(data,trt);
        [p_ics,pout_i]=make_ptable(tmp,'splinefile',options.splinefile,'tablefile',options.tablefile);
        phit_ics=make_phittable(tmp,'splinefile',options.splinefile,'tablefile',options.tablefile);
        disp(['Trait no. ',num2str(trt)]);
        %We need to reshape the 3d tables into 2d matrices, for each
        %trait
        f=@(m) reshape(m,[length(tmp.SNPs.x)*length(cvec),length(logs)]);

        
        trtinf=do_inference(tmp,f(p_ics),f(phit_ics),pout_i);
        inf=concat_inf(inf,trtinf);
    end

else
    disp('SSD and TSD are all we got right now. Choose one.')
end
if options.jack>=0
    pres=get_pres(dataout,inf); 
end
toc;
save(filename,'inf','bootidx','outidx','pres');
end
    
end


function infout=concat_inf(inf1,inf2)

if (~isempty(inf1))&&(~isempty(inf2))
    for fld=fieldnames(inf1)' 
        fld=fld{1};
        infout.(fld)=[inf1.(fld),inf2.(fld)];
    end
elseif isempty(inf1)
    infout=inf2;
elseif isempty(inf2)
    infout=inf1;
else
    infout=[];
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
