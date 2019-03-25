function [MatchedMat,AstOut]= CreateMatchedMatChrom(AstCatalog,varargin)

JD0=2450000;
JD2yr=1/365.25;
Rad2Deg=180/pi;
DefV.ImageSize=[2048 4096];
DefV.ApperaFactor=0.95;
DefV.ImageHighBound=0.8;
DefV.ImageLowBound=0.2;
DefV.MaxMagBound=12;
DefV.MinMagBound=8;
DefV.OnlyAstrometryUsed=false;
DefV.Colls2return={'JD','XWIN_IMAGE','YWIN_IMAGE','MAG_PSF','ALPHAWIN_J2000','DELTAWIN_J2000'};
InPar = InArg.populate_keyval(DefV,varargin,mfilename);

if(InPar.OnlyAstrometryUsed)
    AstCatalog=AstCatOnlyUsedInAstrometry(AstCatalog);
end
for i=1:length(AstCatalog)
    AstCatalog(i).Cat(:,AstCatalog(i).Col.ALPHAWIN_J2000)=AstCatalog(i).Cat(:,AstCatalog(i).Col.ALPHAWIN_J2000)*Rad2Deg;
    AstCatalog(i).Cat(:,AstCatalog(i).Col.DELTAWIN_J2000)=AstCatalog(i).Cat(:,AstCatalog(i).Col.DELTAWIN_J2000)*Rad2Deg;
end

[AstOut,AstUM]=match(AstCatalog);


%Aplly the logical vector to keep onlt the matched objects.
for i=1:length(AstOut)
    %AstOut(i).Cat=AstOut(i).Cat(ConditionVector,:);
    %[RA,Dec]=xy2coo(AstCatalog(i).WCS,[AstOut(i).Cat(:,1) AstOut(i).Cat(:,2)],'OutUnits','deg');
    JD=(cell2mat(AstCatalog(i).getkey('OBSJD')));
    %AstOut(i).Cat=[AstOut(i).Cat RA Dec JD*ones(size(RA))];
    AstOut(i).Cat=[AstOut(i).Cat JD*ones(size(AstOut(i).Cat(:,1)))];
    %AstOut(i).ColCell{end+1}='RA';
    %AstOut(i).ColCell{end+1}='Dec';
    %AstOut(i).ColCell{end+1}='JD';
    %AstOut(i).Col.RA=48;
    %AstOut(i).Col.Dec=49;
    AstOut(i).Col.JD=length(AstOut(i).Cat(1,:));
end
Colls2return=InPar.Colls2return;
[Res,Summary,N_Ep]=astcat2matched_array(AstOut,Colls2return)  ;
ConditionForAppearence=Summary.Nnn>InPar.ApperaFactor*length(AstOut);
ConditionForLocationX= nanmean(Res.XWIN_IMAGE,2)>InPar.ImageLowBound*InPar.ImageSize(1) &nanmean(Res.XWIN_IMAGE,2)<InPar.ImageHighBound*InPar.ImageSize(1);
ConditionForLocationY= nanmean(Res.YWIN_IMAGE,2)>InPar.ImageLowBound*InPar.ImageSize(2) &nanmean(Res.YWIN_IMAGE,2)<InPar.ImageHighBound*InPar.ImageSize(2);
ConditionForMagnitude= nanmean(Res.MAG_PSF,2)<InPar.MaxMagBound &nanmean(Res.MAG_PSF,2)>InPar.MinMagBound;
CondTot=ConditionForAppearence ...
    & ConditionForLocationX & ConditionForLocationY ...
    &ConditionForMagnitude;

for i=1:length(Colls2return)
    MatchedMat.(Colls2return{i})=Res.(Colls2return{i})(CondTot,:);
end

