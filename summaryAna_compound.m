function ana = summaryAna_compound(a_s)
% a_s- the data of all participants

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

% sunnary statistics of main effects
a_sPrev=[nanmean(a_s(:,1:2),2) nanmean(a_s(:,3:4),2)];
a_sCurr=[nanmean(a_s(:,[1 3]),2) nanmean(a_s(:,[2 4]),2)];

mPrev=mean(a_sPrev);
mCurr=mean(a_sCurr);
sePrev=std(a_sPrev)/sqrt(nS); % standard error of the mean
seCurr=std(a_sCurr)/sqrt(nS); % standard error of the mean
ciPrev=[mPrev-1.96*sePrev;mPrev+1.96*sePrev];
ciCurr=[mCurr-1.96*seCurr;mCurr+1.96*seCurr];

% difference scores for the same previous trials
da_s=nan(nS,2);
da_s(:,1)=a_s(:,1)-a_s(:,2);
da_s(:,2)=a_s(:,3)-a_s(:,4);

mda_s=nanmean(da_s);
seda_s=std(da_s)/sqrt(nS);

% difference scores (adaptation effect)
da_sPrev=a_sPrev(:,1)-a_sPrev(:,2);
mda_sPrev=mean(da_sPrev);
seda_sPrev=std(da_sPrev)/sqrt(nS); % standard error of the mean
seda_s_rangePrev=[mda_sPrev-seda_sPrev;mda_sPrev+seda_sPrev];
cida_sPrev=1.96*seda_sPrev; % CI
cida_s_rangePrev=[mda_sPrev-cida_sPrev;mda_sPrev+cida_sPrev];

% difference scores (Pavlovian effect)
da_sCurr=a_sCurr(:,1)-a_sCurr(:,2);
mda_sCurr=mean(da_sCurr);
seda_sCurr=std(da_sCurr)/sqrt(nS); % standard error of the mean
seda_s_rangeCurr=[mda_sCurr-seda_sCurr;mda_sCurr+seda_sCurr];
cida_sCurr=1.96*seda_sCurr; % CI
cida_s_rangeCurr=[mda_sCurr-cida_sCurr;mda_sCurr+cida_sCurr];

% ANOVA
y=a_s;
tbl=table(y(:,1),y(:,2),y(:,3),y(:,4),'VariableNames',{'cAc','sAc','cAs','sAs'});
conds = table(categorical({'c' 'c' 's' 's'}'),categorical({'c' 's' 'c' 's'}'),'VariableNames',{'Prev','Curr'});
rm = fitrm(tbl,'cAc-sAs~1','WithinDesign',conds);
[ranovatbl,~,C,~ ]= ranova(rm,'WithinModel','Prev*Curr');

pPrev=table2array(ranovatbl('(Intercept):Prev','pValue'));
pCurr=table2array(ranovatbl('(Intercept):Curr','pValue'));
pPrevByCurr=table2array(ranovatbl('(Intercept):Prev:Curr','pValue'));

compPrev=multcompare(rm,'Prev');
compCurr=multcompare(rm,'Curr');
compPrevByCurr=multcompare(rm,'Prev','By','Curr','ComparisonType','bonferroni');
compCurrByPrev=multcompare(rm,'Curr','By','Prev','ComparisonType','bonferroni');

ana.indiv=a_s;
ana.m=m;
ana.ci=ci;
ana.se=se;
ana.mPrev=mPrev;
ana.sePrev=sePrev;
ana.ciPrev=ciPrev;
ana.mCurr=mCurr;
ana.seCurr=seCurr;
ana.ciCurr=ciCurr;

ana.da_sPrev=da_sPrev;
ana.mda_sPrev=mda_sPrev;
ana.seda_s_rangePrev=seda_s_rangePrev;
ana.cida_s_rangePrev=cida_s_rangePrev;
ana.da_sCurr=da_sCurr;
ana.mda_sCurr=mda_sCurr;
ana.seda_s_rangeCurr=seda_s_rangeCurr;
ana.cida_s_rangeCurr=cida_s_rangeCurr;

ana.ranovatbl=ranovatbl;
ana.pPrev=pPrev;
ana.pCurr=pCurr;
ana.pPrevByCurr=pPrevByCurr;
ana.compPrevByCurr=compPrevByCurr;
ana.compCurrByPrev=compCurrByPrev;
ana.lilh=lilh;
ana.lilp=lilp;

end

