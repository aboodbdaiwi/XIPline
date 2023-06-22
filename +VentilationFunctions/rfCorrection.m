function NormMR=rfCorrection_djr(MR,maskarray)

Xseg = MR.*(maskarray>0)./(maskarray>0);

for i=1:size(Xseg,1)
    AtP(i)=nanmean(nanmean(Xseg(i,:,:)));
end

for j=1:size(Xseg,2)
    RtL(j)=nanmean(nanmean(Xseg(:,j,:)));
end

for k=1:size(Xseg,3)
    FtH(k)=nanmean(nanmean(Xseg(:,:,k)));
end
corr = Xseg~=Xseg;
corr = transpose(ones(size(Xseg,2),1)*AtP) + ones(size(Xseg,1),1)*RtL;
imshow(corr,[]);

for k=1:size(Xseg,3)
    NormMR(:,:,k)=MR(:,:,k)./(corr*FtH(k));
end
NormMR=NormMR/nanmean(nanmean(nanmean(NormMR)));
NormMR(isnan(NormMR))=0;

return

