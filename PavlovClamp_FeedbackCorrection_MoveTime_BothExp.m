clc
clear
close all

for exp=1:2 % 1- Differential, 2- Compound
    if exp==1
        load('AssociativeAdaptation_Exp1_Differential_trials');
    else
        load('AssociativeAdaptation_Exp2_Compound_trials');
    end
    
    nS=max(T.SN);
    nT=max(T.TN);
    nB=max(T.BN);
    
    clamp_ccw=reshape(T.CCW,nT,[])';
    clamp_ccw_s=clamp_ccw(:,1);
    
    block_mat=reshape(T.BN,nT,[]);
    block=block_mat(:,1);
    
    clamp_i=find(block==1,1);
    wash_i=find(block==2,1);
    
    tr_clamp=clamp_i:(wash_i-1);
    tr_wash=wash_i:nT;
    
    group=reshape(T.group,nT,[])';
    group_s=group(:,1);
    s_test=find(group_s);
    
    % index of 50 ms into the movement
    dt_tablet = 0.005;  % the data is writen in 0.002 sec (500 Hz) and the repeated samples are removed. So the important sampling rate is that of the tablet.
    s2ms=1000;
    i50=50/(dt_tablet*s2ms)+1;
    
    lw=2;
    ms=10;
    
    ha_targ_org=reshape(T.hand_theta,nT,[])';
    ha_50ms_org=reshape(T.hand_theta_50,nT,[])';
    
    ha_targ=ha_targ_org;
    ha_50ms=ha_50ms_org;
    
    ha_targ((clamp_ccw_s==1),:)=-ha_targ((clamp_ccw_s==1),:);
    ha_50ms((clamp_ccw_s==1),:)=-ha_50ms((clamp_ccw_s==1),:);
    
    % Reaction Time
    rt=reshape(T.RT,nT,[])';
    rt_thresh=0.4; % "start faster" message in the experiment
    
    rt_above_thresh_mat=zeros(nS,nT);
    rt_above_thresh_mat(rt>rt_thresh)=1;
    
    tr_high_rt=cell(nS,2); % one column for adaptation and a second for washout
    nTr_high_rt=nan(nS,2); % one column for adaptation and a second for washout
    for s=1:nS
        tr_high_rt{s,1}=find(rt(s,tr_clamp)>rt_thresh); % for each subject
        tr_high_rt{s,2}=find(rt(s,tr_wash)>rt_thresh); % for each subject
        nTr_high_rt(s,1)=length(tr_high_rt{s,1});
        nTr_high_rt(s,2)=length(tr_high_rt{s,2});
    end
    
    % Movement Time
    mt=reshape(T.MT,nT,[])';
    mt_thresh=0.3; % "move faster" message in the experiment
    
    mt_above_thresh_mat=zeros(nS,nT);
    mt_above_thresh_mat(mt>mt_thresh)=1;
    
    tr_high_mt=cell(nS,2); % one column for adaptation and a second for washout
    nTr_high_mt=nan(nS,2); % one column for adaptation and a second for washout
    for s=1:nS
        tr_high_mt{s,1}=find(mt(s,tr_clamp)>mt_thresh); % for each subject
        tr_high_mt{s,2}=find(mt(s,tr_wash)>mt_thresh); % for each subject
        nTr_high_mt(s,1)=length(tr_high_mt{s,1});
        nTr_high_mt(s,2)=length(tr_high_mt{s,2});
    end
    
    ha_targ = removeOutlierTrials(ha_targ,wash_i);
    ha_50ms = removeOutlierTrials(ha_50ms,wash_i);
    
    % remove trials with high reaction time
    ha_targ(find(rt_above_thresh_mat))=nan;
    ha_50ms(find(rt_above_thresh_mat))=nan;
    
    % remove trials with high movement time
    ha_targ(find(mt_above_thresh_mat))=nan;
    ha_50ms(find(mt_above_thresh_mat))=nan;
    
    fb_correct=ha_targ-ha_50ms;
    
    EXP.FBcorr.adapt.all=fb_correct(:,tr_clamp);
    EXP.FBcorr.adapt.m_s=nanmean(EXP.FBcorr.adapt.all,2);
    EXP.FBcorr.adapt.mm_s=nanmean(EXP.FBcorr.adapt.m_s);
    EXP.FBcorr.adapt.sem_s=nanstd(EXP.FBcorr.adapt.m_s)/sqrt(nS);
    EXP.FBcorr.adapt.all_vec=reshape(EXP.FBcorr.adapt.all,[],1);
    EXP.FBcorr.adapt.m_vec=nanmean(EXP.FBcorr.adapt.all_vec);
    EXP.FBcorr.adapt.sd_vec=nanstd(EXP.FBcorr.adapt.all_vec);
    EXP.FBcorr.adapt.se_vec=nanstd(EXP.FBcorr.adapt.all_vec)/sqrt(length(EXP.FBcorr.adapt.all_vec));
    
    mt_afExclude=mt;
    mt_afExclude(find(mt_above_thresh_mat))=nan;
    
    EXP.MT.adapt.all=mt_afExclude(:,tr_clamp);
    EXP.MT.adapt.m_s=nanmean(EXP.MT.adapt.all,2);
    EXP.MT.adapt.mm_s=nanmean(EXP.MT.adapt.m_s);
    EXP.MT.adapt.sem_s=nanstd(EXP.MT.adapt.m_s)/sqrt(nS);
    EXP.MT.adapt.all_vec=reshape(EXP.MT.adapt.all,[],1);
    EXP.MT.adapt.m_vec=nanmean(EXP.MT.adapt.all_vec);
    EXP.MT.adapt.sd_vec=nanstd(EXP.MT.adapt.all_vec);
    EXP.MT.adapt.se_vec=nanstd(EXP.MT.adapt.all_vec)/sqrt(length(EXP.MT.adapt.all_vec));
    
    rt_afExclude=rt;
    rt_afExclude(find(rt_above_thresh_mat))=nan;
    
    EXP.RT.adapt.all=rt_afExclude(:,tr_clamp);
    EXP.RT.adapt.m_s=nanmean(EXP.RT.adapt.all,2);
    EXP.RT.adapt.mm_s=nanmean(EXP.RT.adapt.m_s);
    EXP.RT.adapt.sem_s=nanstd(EXP.RT.adapt.m_s)/sqrt(nS);
    EXP.RT.adapt.all_vec=reshape(EXP.RT.adapt.all,[],1);
    EXP.RT.adapt.m_vec=nanmean(EXP.RT.adapt.all_vec);
    EXP.RT.adapt.sd_vec=nanstd(EXP.RT.adapt.all_vec);
    EXP.RT.adapt.se_vec=nanstd(EXP.RT.adapt.all_vec)/sqrt(length(EXP.RT.adapt.all_vec));
    
    % for the differential exp, separate CS+ and CS- trials
    if exp==1
        % Rotation trials
        rot=reshape(T.ri,nT,[])';

        fb_correct_15clamp=nan(nS,length(tr_clamp)/2);
        fb_correct_0clamp=nan(nS,length(tr_clamp)/2);
        
        for s=1:nS
            fb_correct_15clamp(s,:)=EXP.FBcorr.adapt.all(rot(s,tr_clamp)~=0);
            fb_correct_0clamp(s,:)=EXP.FBcorr.adapt.all(rot(s,tr_clamp)==0);
        end
        
        mfb_correct_15clamp_s=nanmean(fb_correct_15clamp,2);
        mfb_correct_0clamp_s=nanmean(fb_correct_0clamp,2);
        
        fb_correct_15clamp_vec=reshape(fb_correct_15clamp,[],1);
        fb_correct_0clamp_vec=reshape(fb_correct_0clamp,[],1);
        
        mfb_correct_15clamp_vec=nanmean(fb_correct_15clamp_vec);
        mfb_correct_0clamp_vec=nanmean(fb_correct_0clamp_vec);
        
        sdfb_correct_15clamp_vec=nanstd(fb_correct_15clamp_vec);
        sdfb_correct_0clamp_vec=nanstd(fb_correct_0clamp_vec);
        
        sefb_correct_15clamp_vec=nanstd(fb_correct_15clamp_vec)/sqrt(length(fb_correct_15clamp_vec));
        sefb_correct_0clamp_vec=nanstd(fb_correct_0clamp_vec)/sqrt(length(fb_correct_0clamp_vec));
        
        EXP.FBcorr.adapt.fb_correct_15clamp.all=fb_correct_15clamp;
        EXP.FBcorr.adapt.fb_correct_15clamp.m_s=mfb_correct_15clamp_s;
        EXP.FBcorr.adapt.fb_correct_15clamp.mm_s=nanmean(EXP.FBcorr.adapt.fb_correct_15clamp.m_s);
        EXP.FBcorr.adapt.fb_correct_15clamp.sem_s=nanstd(EXP.FBcorr.adapt.fb_correct_15clamp.m_s)/sqrt(nS);
        EXP.FBcorr.adapt.fb_correct_15clamp.all_vec=fb_correct_15clamp_vec;
        EXP.FBcorr.adapt.fb_correct_15clamp.m_vec=mfb_correct_15clamp_vec;
        EXP.FBcorr.adapt.fb_correct_15clamp.sd_vec=sdfb_correct_15clamp_vec;
        EXP.FBcorr.adapt.fb_correct_15clamp.se_vec=sefb_correct_15clamp_vec;
    
        EXP.FBcorr.adapt.fb_correct_0clamp.all=fb_correct_0clamp;
        EXP.FBcorr.adapt.fb_correct_0clamp.m_s=mfb_correct_0clamp_s;
        EXP.FBcorr.adapt.fb_correct_0clamp.mm_s=nanmean(EXP.FBcorr.adapt.fb_correct_0clamp.m_s);
        EXP.FBcorr.adapt.fb_correct_0clamp.sem_s=nanstd(EXP.FBcorr.adapt.fb_correct_0clamp.m_s)/sqrt(nS);
        EXP.FBcorr.adapt.fb_correct_0clamp.all_vec=fb_correct_0clamp_vec;
        EXP.FBcorr.adapt.fb_correct_0clamp.m_vec=mfb_correct_0clamp_vec;
        EXP.FBcorr.adapt.fb_correct_0clamp.sd_vec=sdfb_correct_0clamp_vec;
        EXP.FBcorr.adapt.fb_correct_0clamp.se_vec=sefb_correct_0clamp_vec;
        
        [h,p,ci,stats] = ttest(EXP.FBcorr.adapt.fb_correct_15clamp.m_s,EXP.FBcorr.adapt.fb_correct_0clamp.m_s);
        % Cohen's d effect size
        cohend = computeCohen_d(EXP.FBcorr.adapt.fb_correct_15clamp.m_s, EXP.FBcorr.adapt.fb_correct_0clamp.m_s, 'paired');

        EXP.FBcorr.adapt.ttest_15vs0.h=h;
        EXP.FBcorr.adapt.ttest_15vs0.p=p;
        EXP.FBcorr.adapt.ttest_15vs0.ci=ci;
        EXP.FBcorr.adapt.ttest_15vs0.stats=stats;
        EXP.FBcorr.adapt.ttest_15vs0.cohend=cohend;
        EXP.FBcorr.adapt.ttest_15vs0.mDiff=mean(EXP.FBcorr.adapt.fb_correct_15clamp.m_s-EXP.FBcorr.adapt.fb_correct_0clamp.m_s);
    end
    
    if exp==1
        DIFF=EXP;
    else
        COMP=EXP;
    end
end

FBcorr_BothExp.adapt.all_vec=[DIFF.FBcorr.adapt.all_vec; COMP.FBcorr.adapt.all_vec];
FBcorr_BothExp.adapt.m_vec=nanmean(FBcorr_BothExp.adapt.all_vec);
FBcorr_BothExp.adapt.sd_vec=nanstd(FBcorr_BothExp.adapt.all_vec);
FBcorr_BothExp.adapt.se_vec=nanstd(FBcorr_BothExp.adapt.all_vec)/sqrt(length(FBcorr_BothExp.adapt.all_vec));

MT_BothExp.adapt.all_vec=[DIFF.MT.adapt.all_vec; COMP.MT.adapt.all_vec];
MT_BothExp.adapt.m_vec=nanmean(MT_BothExp.adapt.all_vec);
MT_BothExp.adapt.sd_vec=nanstd(MT_BothExp.adapt.all_vec);
MT_BothExp.adapt.se_vec=nanstd(MT_BothExp.adapt.all_vec)/sqrt(length(MT_BothExp.adapt.all_vec));

RT_BothExp.adapt.all_vec=[DIFF.RT.adapt.all_vec; COMP.RT.adapt.all_vec];
RT_BothExp.adapt.m_vec=nanmean(RT_BothExp.adapt.all_vec);
RT_BothExp.adapt.sd_vec=nanstd(RT_BothExp.adapt.all_vec);
RT_BothExp.adapt.se_vec=nanstd(RT_BothExp.adapt.all_vec)/sqrt(length(RT_BothExp.adapt.all_vec));
    