function C = chrom_comp_stat(GAIAused,GAIAchrom, compused,compchrom)
Rad2miliArcsec= 3600*1000*180/pi;
RA0 = [];
Dec0 =[] ; 
muRA_used = []; 
muDec_used= [];
for i =1:numel(compused.RAfit)
    RA0(i)  =  compused.RAfit(i).fit.p2;
    Dec0(i)    =  compused.Decfit(i).fit.p2;
    muRA_used(i) = compused.RAfit(i).fit.p1*cos(Dec0(i));
    muDec_used(i) = compused.Decfit(i).fit.p1;
end
RA0chrom =  [];
Dec0chrom = [] ; 
muRA_chrom= [];
muDec_chrom=[];
for i =1:numel(compchrom.RAfit)
    RA0chrom(i)  =  compchrom.RAfit(i).fit.p2;
    Dec0chrom(i)    =  compchrom.Decfit(i).fit.p2;
    muRA_chrom(i) = compused.RAfit(i).fit.p1*cos(Dec0chrom(i));
    muDec_chrom(i) = compused.Decfit(i).fit.p1;
end

%{
histogram((RA0'- GAIAused.Cat(:,GAIAused.Col.RA))*Rad2miliArcsec)
title('\alpha original')
xlabel('\Delta  miliarsc');
figure;
histogram((RA0chrom'- GAIAchrom.Cat(:,GAIAchrom.Col.RA))*Rad2miliArcsec)
title('\alpha chromatic')
xlabel('\Delta  miliarsc');
%}
figure; 
plot(Rad2miliArcsec*muRA_used', GAIAused.Cat(:,GAIAused.Col.PMRA),'.');
title('\mu_{\alpha} original vs GAIA')
xlabel('\mu_{\alpha} [marsc]')
ylabel('\mu_{\alpha}^{Gaia} [marsc]')
figure;
plot(Rad2miliArcsec*muRA_chrom', GAIAchrom.Cat(:,GAIAchrom.Col.PMRA),'.');
title('\mu_{\alpha} chrom vs GAIA')
xlabel('\mu_{\alpha} [marsc]')
ylabel('\mu_{\alpha}^{Gaia} [marsc]')

C.RA_std_used = std(RA0'- GAIAused.Cat(:,1))* Rad2miliArcsec;
C.Dec_std_used =std(Dec0'- GAIAused.Cat(:,2))* Rad2miliArcsec;


C.RA_std_chrom = std(RA0chrom'- GAIAchrom.Cat(:,1))* Rad2miliArcsec;
C.Dec_std_chrom =std(Dec0chrom'- GAIAchrom.Cat(:,2))* Rad2miliArcsec;


C.RA_rstd_used = Util.stat.rstd(RA0'- GAIAused.Cat(:,1))* Rad2miliArcsec;
C.Dec_rstd_used =Util.stat.rstd(Dec0'- GAIAused.Cat(:,2))* Rad2miliArcsec;

C.RA_rstd_chrom = Util.stat.rstd(RA0chrom'- GAIAchrom.Cat(:,1))* Rad2miliArcsec;
C.Dec_rstd_chrom =Util.stat.rstd(Dec0chrom'- GAIAchrom.Cat(:,2))* Rad2miliArcsec;

C.mu_RA_rstd_used = Util.stat.rstd(Rad2miliArcsec*muRA_used'- GAIAused.Cat(:,GAIAused.Col.PMRA));
C.mu_Dec_rstd_used = Util.stat.rstd(Rad2miliArcsec*muDec_used'- GAIAused.Cat(:,GAIAused.Col.PMDec));

C.mu_RA_rstd_chrom = Util.stat.rstd(Rad2miliArcsec*muRA_chrom'- GAIAchrom.Cat(:,GAIAchrom.Col.PMRA));
C.mu_Dec_rstd_chrom = Util.stat.rstd(Rad2miliArcsec*muDec_chrom'- GAIAchrom.Cat(:,GAIAchrom.Col.PMDec));

end