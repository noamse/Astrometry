PTFDirectory='/data/noamse/Astrometry/Data/Catalogs/catPTF_6_1_19/';
ZTFDirectory='/data/noamse/Astrometry/Data/Catalogs/ZTF_14_3_19/';

[PTF_GAIAchrom,PTF_GAIAused,PTF_datachrom, PTF_dataused, PTF_compchrom,PTF_compused,PTF_astcatused,PTF_Cat]= chromatic_comp(PTFDirectory);
%[ZTF_GAIAchrom,ZTF_GAIAused,ZTF_datachrom, ZTF_dataused, ZTF_compchrom,ZTF_compused,ZTF_astcatused,ZTF_Cat]= chromatic_comp(ZTFDirectory,'Survey','ZTF');
PTFstat = chrom_comp_stat(PTF_GAIAused,PTF_GAIAchrom, PTF_compused,PTF_compchrom);
%ZTFstat = chrom_comp_stat(ZTF_GAIAused,ZTF_GAIAchrom, ZTF_compused,ZTF_compchrom);
%{
%%

PTFDirectory='/data/noamse/Astrometry/Data/Catalogs/catPTF_7_3_19/';
astcatused = Asmtry.open_directory_astcat('Directory',PTFDirectory);
astcatused = clear_failure(astcatused);
astcatused = Asmtry.get_astrometry_res(astcatused);
[CatZTF,astcatchrom]= pa_chromatic_corr(astcatused);
[dataused,~]= Asmtry.match_mat(astcatused,'match_SearchRadius',2);
[datachrom,~]= Asmtry.match_mat(astcatchrom,'match_SearchRadius',2);
[GAIAused, compused]   = Asmtry.compare_astrometry_gaia(astcatused,dataused);
[GAIAchrom, compchrom] = Asmtry.compare_astrometry_gaia(astcatchrom,datachrom);



%%
RA0used = [];
Dec0used =[] ; 
for i =1:numel(compused.RAfit)
    RA0used(i)  =  compused.RAfit(i).fit.p2;
    Dec0used(i)    =  compused.Decfit(i).fit.p2;
end
RA0chrom = [];
Dec0chrom =[] ; 
for i =1:numel(compchrom.RAfit)
    RA0chrom(i)  =  compchrom.RAfit(i).fit.p2;
    Dec0chrom(i)    =  compchrom.Decfit(i).fit.p2;
end


histogram(RA0used'- GAIAused.Cat(:,1))
histogram(RA0chrom'- GAIAchrom.Cat(:,1))

C.RAstd_used = std(RA0used'- GAIAused.Cat(:,1))* 3600*180/pi*1000;
C.Decstd_used =std(Dec0used'- GAIAused.Cat(:,2))* 3600*180/pi*1000;

C.RAstd_chrom = std(RA0chrom'- GAIAchrom.Cat(:,1))* 3600*180/pi*1000;
C.Decstd_chrom =std(Dec0chrom'- GAIAchrom.Cat(:,2))* 3600*180/pi*1000;
%}