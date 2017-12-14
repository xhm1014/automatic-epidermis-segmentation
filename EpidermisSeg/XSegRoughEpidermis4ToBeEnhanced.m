%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this is a function used to find  the partly rough epidermis area,
% for the enhanced method. Entact epidermis is not required, but highly 
% highly true positive rate is required.

% Input:
%   -IM    RGB image
%   -Channel  indicates which channel will be used
% Output:
%   -maskEpidermis    a logical matrix indicate the position of the
%                   epidermis
% Key Threshold:
%   -TAxisRatio % we define the enlonged one as the AxisRatio > TAxisRatio
%   -TAreaofROI the threshold for the size of area that we think is noise

% (c) Edited by Cheng Lu,
% Deptment of Eletrical and Computer Engineering,
% University of Alberta, Canada.  20th Feb, 2010
% If you have any problem feel free to contact me.
% Please address questions or comments to: hacylu@yahoo.com
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function maskEpidermis=LSegRoughEpidermis4ToBeEnhanced(IM,Channel,shown)

sizeIM=size(IM);
% predefined area threshold
TAreaofROI=ceil(sizeIM(1)*sizeIM(2)/1000);
%TAreaofROI=ceil(sizeIM(1)*sizeIM(2)/500);
% show(IM);
%  in R,G,B channels, R is the better one, so
%%  Thresholding in predefined channel
if ~exist('Channel','var')||strcmp(Channel,'R')
    Rchannel=IM(:,:,1);
else
    
    if strcmp(Channel,'G')
        Rchannel=IM(:,:,2);
    end
    
    if strcmp(Channel,'B')
        Rchannel=IM(:,:,3);
    end
    
    if strcmp(Channel,'Gray')
        Rchannel=rgb2gray(IM);
    end
    
    if strcmp(Channel,'H')
        HSV=rgb2hsv(IM);
        Rchannel=HSV(:,:,1);
    end
    
    if strcmp(Channel,'S')
        HSV=rgb2hsv(IM);
        Rchannel=HSV(:,:,2);
    end
    
    if strcmp(Channel,'V')
        HSV=rgb2hsv(IM);
        Rchannel=HSV(:,:,3);
    end
end

if ~exist('shown','var')
    shown=0;
end
% show(Rchannel);
% Rchannel=adapthisteq(Rchannel); %% added for CLAHE enhancing image contrast
%% segment out the background first
TBKWhite=0.95;
RCnoBK=im2bw(Rchannel,TBKWhite);
RCnoBK=~RCnoBK;

% show(RCnoBK,11);
%% thresholding in the foreground
IM_Temp=IM(RCnoBK);
% hist(double(IM_Temp),255);
T4IM=graythresh(IM_Temp);% T4IM*255
% T4IM=T4IM; % in order to focus on the epidermis
IMlogical= im2bw(Rchannel,T4IM);
%         IMlogical= IM<T4IM;%
RC_Thresh1=~IMlogical;
if shown
    show(RC_Thresh1);
end
% RC_Thresh1=LThresholding(Rchannel);

% good for M1,M4,M115
% RC_Thresh1=LRecursiveThresholding(Rchannel,2);

% good for M1 M4 but may slow the process
% RC_Thresh1=LThresholding(Rchannel);
% show(RC_Thresh1,12)
% show(IM)
%% 
% SE=strel('disk',2);
% RC_Thresh1=imerode(RC_Thresh1,SE);
% % show(RC_Thresh1)
RC_Thresh1_open=bwareaopen(RC_Thresh1, TAreaofROI,4);
if shown
    show(RC_Thresh1_open);
end
%% no dilate
% % RC_Thresh1_dilate=imdilate(RC_Thresh1,SE);
% % show(RC_Thresh1_dilate)
% RC_Thresh1_open=bwareaopen(RC_Thresh1, TAreaofROI,4);
% % show(RC_Thresh1_open,13);
%
% % tic
% % se=strel('disk',5,8);
% % RC_Thresh1_open_v2=imopen(RC_Thresh1,se);
% % toc
%
% % show(RC_Thresh1_open_v2);
% % show(RC_Thresh1_open);
%%
% measuring the big component(assume only 1)
% CC=bwconncomp(RC_Thresh1_open);
% [L,numObjs]=bwlabel(RC_Thresh1_open);


%% analysis the remaining objs, find the longest one
CC=bwconncomp(RC_Thresh1_open);
% TAxisRatio=4; % we define the enlonged one as the AxisRatio >3 for special case
TAxisRatio=2; % new test parameter
idx4Candidate=1;
CandidateSet=[];
if CC.NumObjects>1
    STATStemp=regionprops(CC,'MajorAxisLength','MinorAxisLength','Perimeter');
    AxisRatiotemp=[STATStemp.MajorAxisLength]./[STATStemp.MinorAxisLength];
    Perimetertemp=[STATStemp.Perimeter];
       
%     TAvePremeter=mean(Perimetertemp);
    
     Max_Premeter=max(Perimetertemp);
    for i=1:CC.NumObjects
        if AxisRatiotemp(i)>TAxisRatio && Perimetertemp(i)>(1/3)*Max_Premeter
%        if AxisRatiotemp(i)<4  | AxisRatiotemp(i)>10                                       %%%%%%%%%%% for special test case
            CandidateSet(idx4Candidate)=i;
            idx4Candidate=idx4Candidate+1;
        end
    end
    
    if isempty(CandidateSet)
        % new weighted judgement
        AxisRatiotempN=LNorHist(AxisRatiotemp);
        PerimetertempN=LNorHist(Perimetertemp);
        TempPlus=AxisRatiotempN+PerimetertempN;
        [AxisRatiotemp_max,CandidateSet]=max(TempPlus);
        % original code
%         [AxisRatiotemp_max,CandidateSet]=max(AxisRatiotemp);
        List_RegionwithMaxAxisRatio=CC.PixelIdxList{CandidateSet};
    else
        List_RegionwithMaxAxisRatio=[];
        for i=1:length(CandidateSet)
            List_RegionwithMaxAxisRatio=[List_RegionwithMaxAxisRatio;CC.PixelIdxList{CandidateSet(i)}];
        end
    end
end

if CC.NumObjects==1
        List_RegionwithMaxAxisRatio=CC.PixelIdxList{1};
end

if CC.NumObjects==0
        error('There is no Epidermis?????');
end

clear STATStemp MaxAxisRatio CandidateSet;
maskEpidermis=RC_Thresh1_open;
maskEpidermis(:)=0;
maskEpidermis(List_RegionwithMaxAxisRatio)=1;

end