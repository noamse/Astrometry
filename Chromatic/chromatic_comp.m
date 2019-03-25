function [GAIAchrom,GAIAused,datachrom, dataused, compchrom,compused,astcatused,Cat]= chromatic_comp(Directory,varargin)

DefV.Survey ='PTF';
DefV.SearchRadius= 1;
InPar = InArg.populate_keyval(DefV,varargin,mfilename);


astcatused = Asmtry.open_directory_astcat('Directory',Directory);
astcatused = clear_failure(astcatused);
%astcatused = Asmtry.units_change(astcatused);
astcatused = Asmtry.get_astrometry_res(astcatused);

[Cat,astcatchrom]= pa_chromatic_corr(astcatused,'Survey',InPar.Survey);
%[Cat,astcatchrom]= pa_chromatic_corr(astcatused);

[dataused,~]= Asmtry.match_mat(astcatused,'match_SearchRadius',2,'Survey',InPar.Survey);
[datachrom,~]= Asmtry.match_mat(astcatchrom,'match_SearchRadius',2,'Survey',InPar.Survey);


[GAIAused, compused]   = Asmtry.compare_astrometry_gaia(astcatused,dataused);
[GAIAchrom, compchrom] = Asmtry.compare_astrometry_gaia(astcatchrom,datachrom);
end
