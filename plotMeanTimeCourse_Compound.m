function plotMeanTimeCourse_Compound(a_s,tr_b)
% a- the data of all participants
% g- group to analyze (1- conditioning; 0- control in g_s)
% s- subjects to analyze

set(groot, 'DefaultAxesXColor', [0,0,0], ...
           'DefaultAxesYColor', [0,0,0], ...
           'DefaultAxesZColor', [0,0,0])
       
co= [174,89,116]/255; % pink

nS=size(a_s,1); % number of subjects;
nT=size(a_s,2); % number of trials;
wash_i=tr_b(1);

% bootstrap for non-parametric confidence interval of the mean
nboot=1000;
bootfun=@(x)nanmean(x);

% mean and ci of each condition (each cell)
ma_s=nanmean(a_s);
cia_s=bootci(nboot,bootfun,a_s);
sea_s=nanstd(a_s)/sqrt(nS);

% vara_s=cia_s;
vara_s=[ma_s-sea_s;ma_s+sea_s];

c_noFb=[1 1 1]*0.9;

lw=5;
afs=26;
tfs=22;

xLimC=[0 nT];
yLim_HA=[-5 35];

xTick=0:100:900;
yTick_HA=-80:10:80;

txt_yloc=yLim_HA(2)-5;

trA=1:(wash_i-1);
trB=wash_i:nT;

figure('position',[50 100 2025 400])
hold on
patch('Faces',1:4,'Vertices',[wash_i yLim_HA(1); nT yLim_HA(1); nT+1 yLim_HA(2); wash_i yLim_HA(2)],'FaceColor',c_noFb,'EdgeColor','none')
plot(ones(1,2)*(wash_i-.5),yLim_HA,':','color',[1 1 1]*0.3,'linewidth',3)
plot(xLimC,[0 0],':k','linewidth',2)
fill([trA flip(trA)],[vara_s(1,trA) flip(vara_s(2,trA))]',co,'linestyle','none','facealpha',0.3);
fill([trB flip(trB)],[vara_s(1,trB) flip(vara_s(2,trB))]',co,'linestyle','none','facealpha',0.3);
plot(trA,ma_s(trA),'color',co,'linewidth',lw)
plot(trB,ma_s(trB),'color',co,'linewidth',lw)
xlabel('Trial Number','fontsize',afs)
% ylabel('Hand Angle (deg)','fontsize',afs)
ylabel('Heading Angle (deg)','fontsize',afs)
set(gca,'xtick',xTick,'ytick',yTick_HA,'xlim',xLimC,'ylim',yLim_HA,'fontsize',tfs,'tickdir','out')



end

