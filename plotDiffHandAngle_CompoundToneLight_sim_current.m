function plotDiffHandAngle_CompoundToneLight_sim_current(a)
% a- the data of all participants
% g- group to analyze (1- conditioning; 0- control in g_s)
% s- subjects to analyze

set(groot, 'DefaultAxesXColor', [0,0,0], ...
    'DefaultAxesYColor', [0,0,0], ...
    'DefaultAxesZColor', [0,0,0])

co_afHit= [174,89,116]/255; % moderate pink


a_s=a;

nS=size(a_s,1); % number of participants;
nC=size(a_s,2); % number of conditions;

% mean and ci of each condition (each cell)
ma_s=nanmean(a_s);

bw=.5;
ms=12;
blocx=1:3;

lw=2;
afs=24;
tfs=18;

xLim=[.5 3.5];

yTick=-6:.5:6;
% yLabel='\DeltaHand Angle (deg)';
yLabel='\DeltaHeading Angle (deg)';

yBase=mean(ma_s);
yLim=[-0.6 0.7];
% yLim=[yBase-0.6 yBase+0.7];

figure('position',[50 100 274 400])
hold on
plot(xLim,[0 0],':','color',[1 1 1]*.2,'linewidth',2)
% plot(xLim,yBase*[1 1],':','color',[1 1 1]*.2,'linewidth',2)

for c=1:nC
    co=co_afHit;
    
    if c==1
        bar(blocx(c),ma_s(c),'basevalue',0,'barwidth',bw,'linewidth',lw, 'edgecolor',co,'facecolor',co);
    else
        bar(blocx(c),ma_s(c),'basevalue',0,'barwidth',bw,'linewidth',lw, 'edgecolor',co,'facecolor','none');
    end
    
%     if c==1
%         bar(blocx(c),ma_s(c),'basevalue',yBase,'barwidth',bw,'linewidth',lw, 'edgecolor',co,'facecolor',co);
%     else
%         bar(blocx(c),ma_s(c),'basevalue',yBase,'barwidth',bw,'linewidth',lw, 'edgecolor',co,'facecolor','none');
%     end
    
end

set(gca,'xtick',blocx,'xticklabel',{'Comp','Tone','Light'},'ytick',yTick,'fontsize',tfs,'tickdir','out')
% set(gca,'xtick',blocx,'xticklabel',{'Comp','Tone','Light'},'ytick',[],'fontsize',tfs,'tickdir','out')
xtickangle(-60)
ylabel(yLabel,'fontsize',afs)
xlim(xLim);
ylim(yLim);

end

