function [] = OpenImageAndAstrometry(varargin)
%call for mextractor and astromery
%Save AstCat files with some options.
DefV.UseCase_TranC      = {'affine_tt_cheby2_4', 100; 'affine_tt_cheby2_3', 70; 'affine_tt',          10; 'affine',             5};
DefV.SaveDirectory='/home/noamse/matlab/rms_measure/OpenedImages/';
DefV.SaveNameAstCat='AstCat_1.mat';
DefV.FileSize=1;
DefV.StartIndex=1;
DefV.SourceDirectory='/home/noamse/matlab/rms_measure/TestIm/';
DefV.MaxPMerr=[];
DefV.SubImSize=[1500 1500];
DefV.Subimage=false;
DefV.ApplyParallax= true;
DefV.Index=  [];
DefV.Survey= 'PTF';
InPar = InArg.populate_keyval(DefV,varargin,mfilename);

A= dir([InPar.SourceDirectory '*.fits']); 

FileNames=cell(length(A),1);
if isempty(InPar.Index)
   Index = 1:numel(FileNames );
else
    Index = InPar.Index;
end
FileNames={A(Index).name};  

for i=1:numel(FileNames)
    FileNames{i}=strcat(InPar.SourceDirectory,FileNames{i});
end

for i=1:numel(FileNames)
    i
    FileNames{i}
    S = FITS.read2sim(FileNames{i}); 
    S = mextractor(S); 
    switch InPar.Survey
        case 'PTF'
            [R,Sa] = astrometry(S,'UseCase_TranC',InPar.UseCase_TranC,'MaxPMerr',InPar.MaxPMerr,'MaxExcessNoise',1,'ApplyParallax',true); 
        case 'ZTF'
            [R,Sa] = astrometry(S,'UseCase_TranC',InPar.UseCase_TranC,'MaxPMerr',InPar.MaxPMerr,'MaxExcessNoise',1,'ApplyParallax',true,'SCALE',1.012,'Flip',[1 1]); 
    end
    %updating the coordinates in the catalog
    Sa=update_coordinates(Sa);
    %creating a AstCat object to save more efficiently
    AstCatTemp=AstCat.sim2astcat(Sa);
    AstCatTemp.UserData.R=R;
    AstCatTemp.UserData.FileName=FileNames{i};
    
    save([InPar.SaveDirectory InPar.SaveNameAstCat],'AstCatTemp');
    InPar.SaveNameAstCat=['AstCat_' num2str(i) '.mat'];
    
end
end

