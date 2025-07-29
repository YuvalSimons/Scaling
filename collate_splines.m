function collate_splines(spline_folder,splinefile)
logs=(-6:1/16:-1)'; 
for i=1:length(logs)
    filename=[spline_folder,'/',num2str(i),'.mat'];
    S=load(filename);
    sp(i)=S.sp;
    spv(i)=S.spv;
end
save(splinefile,'sp','spv');
end