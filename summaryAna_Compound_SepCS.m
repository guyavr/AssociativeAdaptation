function ana = summaryAna_Compound_SepCS(a_s)
% a- the data of all participants
% s- participants to analyze

nS=size(a_s,1);
nC=size(a_s,2);

% bootstrap for non-parametric confidence interval of the mean
nboot=1000;
bootfun=@(x)nanmean(x);

lilh=nan(1,nC);
lilp=nan(1,nC);

if nS>=4
    for c=1:nC
        [lilh(c),lilp(c)] = lillietest(a_s(:,c));
    end
end
            
m=nanmean(a_s);
ci=bootci(nboot,bootfun,a_s);
se=std(a_s)/sqrt(nS); % standard error of the mean
            
% ANOVA
y=a_s;
tbl=table(y(:,1),y(:,2),y(:,3),'VariableNames',{'comp','tone','frame'});
conds = table(categorical({'c' 't' 'f'}'),'VariableNames',{'Cond'});
rm = fitrm(tbl,'comp-frame~1','WithinDesign',conds);
[ranovatbl,~,C,~ ]= ranova(rm,'WithinModel','Cond');

pCond=table2array(ranovatbl('(Intercept):Cond','pValue'));

mltcmp=multcompare(rm,'Cond','ComparisonType','bonferroni');
% mltcmp=multcompare(rm,'Cond','ComparisonType','tukey-kramer');
% mltcmp=multcompare(rm,'Cond','ComparisonType','lsd');

% correlation
a_tl=a_s(:,2:3);
[rho_pear,pval_pear] = corr(a_tl);
corrPear_r = rho_pear(1,2);
corrPear_p = pval_pear(1,2);

[rho_spear,pval_spear] = corr(a_tl,'Type','Spearman');
corrSpear_r = rho_spear(1,2);
corrSpear_p = pval_spear(1,2);

ana.indiv=a_s;
ana.m=m;
ana.ci=ci;
ana.se=se;
ana.ranovatbl=ranovatbl;
ana.pCond=pCond;
ana.mltcmp=mltcmp;
ana.lilh=lilh;
ana.lilp=lilp;

ana.corrPear_r=corrPear_r;
ana.corrPear_p=corrPear_p;
ana.corrSpear_r=corrSpear_r;
ana.corrSpear_p=corrSpear_p;

end

