function plotBarsMeanDiffHandAngle_CompoundToneLight(a,tr_c)
% a- the data of all participants
% s- subjects to analyze

set(groot, 'DefaultAxesXColor', [0,0,0], ...
           'DefaultAxesYColor', [0,0,0], ...
           'DefaultAxesZColor', [0,0,0])
       
co_afHit= [174,89,116]/255; % moderate pink

a_s=a.indiv;

nS=size(a_s,1); % number of participants;
nC=size(tr_c,2); % number of conditions;

% bootstrap for non-parametric confidence interval of the mean
nboot=1000;
bootfun=@(x)nanmean(x);

% mean and ci of each condition (each cell)
ma_s=nanmean(a_s);
cia_s=bootci(nboot,bootfun,a_s);
sea_s=std(a_s)/sqrt(nS);

% vara_s=cia_s;
vara_s=[ma_s-sea_s;ma_s+sea_s];

bw=.5;
ms=12;
blocx=1:3;

lw=2;
afs=24;
tfs=18;

xLim=[.5 3.5];
yLim=[-0.6 0.7];

yTick=-6:.5:6;
yLabel='\DeltaHeading Angle (deg)';

figure('position',[50 100 300 400])
hold on
plot(xLim,[0 0],':','color',[1 1 1]*.2,'linewidth',2)
for c=1:nC
    co=co_afHit;
    
    if c==1
        bar(blocx(c),ma_s(c),'basevalue',0,'barwidth',bw,'linewidth',lw, 'edgecolor',co,'facecolor',co);
    else
        bar(blocx(c),ma_s(c),'basevalue',0,'barwidth',bw,'linewidth',lw, 'edgecolor',co,'facecolor','none');
    end
    
    plot([blocx(c) blocx(c)],vara_s(:,c),'-k','linewidth',2)
    
end

% add p values
txtfs=16;
xtxtloc=xLim(1)+.5*(xLim(2)-xLim(1));
ydisttxt=.1*(yLim(2)-yLim(1));
ydistline=.06*(yLim(2)-yLim(1));
ytxt_loc1=yLim(2)-0.05;

if a.mltcmp.pValue(1)<0.05 %compound vs frame
    if a.mltcmp.pValue(1)<0.001
        text(xtxtloc,ytxt_loc1,'p<0.001','fontsize',txtfs,'horizontalalignment','center','fontweight','bold')
    else
        text(xtxtloc,ytxt_loc1,sprintf('p=%.3f',a.mltcmp.pValue(1)),'fontsize',txtfs,'horizontalalignment','center','fontweight','bold')
    end
else
    text(xtxtloc,ytxt_loc1,sprintf('p=%.3f',a.mltcmp.pValue(1)),'fontsize',txtfs,'horizontalalignment','center')
end
plot([1 3],ones(1,2)*ytxt_loc1-ydistline,'-k','linewidth',3)

if a.mltcmp.pValue(2)<0.05 %compound vs frame
    if a.mltcmp.pValue(2)<0.001
        text(xtxtloc,ytxt_loc1,'p<0.001','fontsize',txtfs,'horizontalalignment','center','fontweight','bold')
    else
        text(1.5,ytxt_loc1-1.3*ydisttxt,sprintf('p=%.3f',a.mltcmp.pValue(2)),'fontsize',txtfs,'horizontalalignment','center','fontweight','bold')
    end
else
    text(xtxtloc,ytxt_loc1-1.3*ydisttxt,sprintf('p=%.3f',a.mltcmp.pValue(2)),'fontsize',txtfs,'horizontalalignment','center')
end
plot([1 2],ones(1,2)*ytxt_loc1-1.3*ydisttxt-ydistline,'-k','linewidth',3)

set(gca,'xtick',blocx,'xticklabel',{'Comp','Tone','Light'},'ytick',yTick,'fontsize',tfs,'tickdir','out')
xtickangle(-60)
ylabel(yLabel,'fontsize',afs)
xlim(xLim);
ylim(yLim);

a_tl=a_s(:,2:3);
% Plot correlation
xLim_corr=[-1.5 1];
yLim_corr=xLim_corr;
[rho,pval] = corr(a_tl);
corr_r = rho(1,2);
corr_p = pval(1,2);
figure('position',[50 100 400 400])
hold on
plot([-1.5 1],[1 -1.5],':','color',[1 1 1]*.4,'linewidth',3)
pc=scatter(a_tl(:,1),a_tl(:,2),'o','markeredgecolor','none','markerfacecolor',co_afHit,'sizedata',200);
alpha 0.7
uistack(pc,'top')
xlabel({'\DeltaHeading Angle (deg)','Tone Alone'},'fontsize',20)
ylabel({'\DeltaHeading Angle (deg)','Light Alone'},'fontsize',20)
set(gca,'xlim',xLim_corr,'ylim',yLim_corr,'xtick',-1:1:1,'ytick',-1:1:1,'fontsize',18,'tickdir','out')

end

