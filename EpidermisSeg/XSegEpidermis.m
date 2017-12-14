%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this is a function used to segment the epidermis.

% Input:
%   -IM    RGB image
%   -epi_thick epidermis thickness
%   -initial_epidermis_mask  the initially segmented epidermis region
% Output:
%   -maskEpidermis    a logical matrix indicate the position of the
%                   epidermis
% Key Threshold:
%   -TTickness
%   -TAxisRatio % we define the enlonged one as the AxisRatio > TAxisRatio
%   -TAreaofROI the threshold for the size of area that we think is noise

% (c) Edited by Hongming Xu,
% Deptment of Eletrical and Computer Engineering,
% University of Alberta, Canada.  20th Feb, 2010
% If you have any problem feel free to contact me.
% Please address questions or comments to: mxu@ualberat.ca
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [epidermis_mask imagethick]=XSegEpidermis(IM,epi_thick,initial_epidermis_mask,shown)
if ~exist('shown','var')
    shown=0;
end
TTickness=150;
% TTickness=160*(size(IM,1)/1800); % 150 is empirically selected threshold
if epi_thick < TTickness 
    epidermis_mask=initial_epidermis_mask;
    imagethick=epi_thick;
else
    
%      lab = makecform('srgb2lab');
%      IM = applycform(IM,lab);
%      a1=double(IM(:,:,2));b1=double(IM(:,:,3));
%      ind=find(initial_epidermis_mask);
%      a=a1(ind);b=b1(ind);
%      ab=[a,b];
%      [IDX,C]=kmeans(ab,2);
%    %% Gaussian mixture model
%    ind0=find(initial_epidermis_mask);
%    R0=IM(:,:,1);G0=IM(:,:,2);B0=IM(:,:,3);
%    R1=R0(ind0);
%    G1=G0(ind0);
%    B1=B0(ind0);
%    D=double([R1';G1';B1']);
%    nbStates=2;
%    [Priors, Mu, Sigma] = EM_init_kmeans(D, nbStates);
%    [Priors, Mu, Sigma,Pix] = EM(D, Priors, Mu, Sigma);
%    [C,Index]=max(Pix,[],2);
%     ind=find(Index==1);
%     mask1=zeros(size(R0));
%     mask1(ind0(ind))=1;
%     mask2=bitxor(mask1,initial_epidermis_mask);
%     
%     sizeIM=size(initial_epidermis_mask);    
%     TAreaofROI=ceil(sizeIM(1)*sizeIM(2)/1000);% predefined area threshol
%     epidermis_mask1=bwareaopen(mask1,TAreaofROI,4);
%     if shown
%         show(epidermis_mask1,2);
%     end
%     epidermis_mask2=bwareaopen(mask2,TAreaofROI,4);
%     if shown
%         show(epidermis_mask2,3);
%     end

    sizeIM=size(initial_epidermis_mask);    
    TAreaofROI=ceil(sizeIM(1)*sizeIM(2)/1000);% predefined area threshold
%    TAreaofROI=ceil(sizeIM(1)*sizeIM(2)/5000);  %%%%%%%%%%% predefined area threshold for test special case
%% K-means algorithm
    red=IM(:,:,1);green=IM(:,:,2);blue=IM(:,:,3);
    ind=find(initial_epidermis_mask);
    r=red(ind);g=green(ind);b=blue(ind);
    X=double([r,g,b]);
    [IDX,C]=kmeans(X,2,'Replicates',1);   % kmeans clustering algorithm
    
    mask=zeros(size(initial_epidermis_mask));
    if sum(C(1,:))<sum(C(2,:)) %% based on the knowledge epidermis is darker
        id1=find(IDX==1);
        mask(ind(id1))=1;
    else
        id2=find(IDX==2);
        mask(ind(id2))=1;
    end
    
                                 
    RC_Thresh1_open1=bwareaopen(mask, TAreaofROI,8);   
    
    if shown
        ind1=find(RC_Thresh1_open1==0);
        red(ind1)=255;green(ind1)=255;blue(ind1)=255;
        img=cat(3,red,green,blue);
        show(img);
    end
%%--K-mean--%%

%     IM=adapthisteq(IM);  %% added for CLAHE enhancing image contrast

%%-- global thresholding by using Rchannel --%%
%     IM_R=(1/3)*(IM(:,:,1)+IM(:,:,2)+IM(:,:,3));
%     IM_Temp=IM_R(initial_epidermis_mask);
%     ThIM=graythresh(IM_Temp);
%     IMlogical=im2bw(IM_R,ThIM);
%     RC_Thresh1_open1=bwareaopen(~IMlogical, TAreaofROI,8);

    
    %%-- global thresholding by using RGB channels --%%
%     IM_R=IM(:,:,1); IM_G=IM(:,:,2); IM_B=IM(:,:,3);
%     IM_Temp(:,:,1)=IM_R(initial_epidermis_mask);
%     IM_Temp(:,:,2)=IM_G(initial_epidermis_mask);
%     IM_Temp(:,:,3)=IM_B(initial_epidermis_mask);
%     ThIM=multithresh(IM_Temp);
%     imrgb=imquantize(IM,ThIM,[0 255]);
%     Temp=imrgb(:,:,1)+imrgb(:,:,2)+imrgb(:,:,3);
%     ind=find(Temp==0);
%     IMlogical=zeros(size(Temp));
%     IMlogical(ind)=1;
%     RC_Thresh1_open1=bwareaopen(IMlogical, TAreaofROI,8);
%%-- filter the false epidermis regions --%%

     epidermis_mask=zeros(size(initial_epidermis_mask));
     CC=bwconncomp(RC_Thresh1_open1);
     TAxisRatio=3;  % we define the enlonged one as the AxisRatio >3
   
     if CC.NumObjects>1
         STATStemp=regionprops(CC,'MajorAxisLength','MinorAxisLength');
         AxisRatiotemp=[STATStemp.MajorAxisLength]./[STATStemp.MinorAxisLength];
         for i=1:CC.NumObjects
             if AxisRatiotemp(i)>TAxisRatio
                 epidermis_mask(CC.PixelIdxList{i})=1;
             end
         end
         if sum(epidermis_mask(:))==0
            [AxisRatiotemp_max,CandidateSet]=max(AxisRatiotemp);
             epidermis_mask(CC.PixelIdxList{CandidateSet})=1;
         end
     else
         epidermis_mask=RC_Thresh1_open1;
     end
     clear STATStemp TAxisRatio IM initial_epidermis_mask

     imagethick=XThicknessCal(epidermis_mask);
     
     
%      %% analysis the remaining objs, find the longest one
%     epidermis_mask=zeros(size(initial_epidermis_mask));
%      CC=bwconncomp(RC_Thresh1_open1);
%      % TAxisRatio=4; % we define the enlonged one as the AxisRatio >3
%      TAxisRatio=2; % new test parameter
%      idx4Candidate=1;
%      CandidateSet=[];
%      if CC.NumObjects>1
%          STATStemp=regionprops(CC,'MajorAxisLength','MinorAxisLength','Perimeter');
%          AxisRatiotemp=[STATStemp.MajorAxisLength]./[STATStemp.MinorAxisLength];
%          Perimetertemp=[STATStemp.Perimeter];
%          
%          %     TAvePremeter=mean(Perimetertemp);
%          
%          Max_Premeter=max(Perimetertemp);
%          for i=1:CC.NumObjects
%              if AxisRatiotemp(i)>TAxisRatio && Perimetertemp(i)>(0.6)*Max_Premeter
%                  CandidateSet(idx4Candidate)=i;
%                  idx4Candidate=idx4Candidate+1;
%              end
%          end
%          
%          if isempty(CandidateSet)
%              % new weighted judgement
%              AxisRatiotempN=LNorHist(AxisRatiotemp);
%              PerimetertempN=LNorHist(Perimetertemp);
%              TempPlus=AxisRatiotempN+PerimetertempN;
%              [AxisRatiotemp_max,CandidateSet]=max(TempPlus);
%              % original code
%              %         [AxisRatiotemp_max,CandidateSet]=max(AxisRatiotemp);
%              List_RegionwithMaxAxisRatio=CC.PixelIdxList{CandidateSet};
%          else
%              List_RegionwithMaxAxisRatio=[];
%              for i=1:length(CandidateSet)
%                  List_RegionwithMaxAxisRatio=[List_RegionwithMaxAxisRatio;CC.PixelIdxList{CandidateSet(i)}];
%              end
%          end
%      end
%      
%      if CC.NumObjects==1
%          List_RegionwithMaxAxisRatio=CC.PixelIdxList{1};
%      end
%      
%      if CC.NumObjects==0
%          error('There is no Epidermis?????');
%      end
%      
%      clear STATStemp MaxAxisRatio CandidateSet;
%      
%      epidermis_mask(List_RegionwithMaxAxisRatio)=1;
%     
%     if CC.NumObjects>1
%         epidermis_mask1=logical(epidermis_mask);
%         red=immultiply(epidermis_mask1,IM(:,:,1));
%         green=immultiply(epidermis_mask1,IM(:,:,2));
%         blue=immultiply(epidermis_mask1,IM(:,:,3));
%         g=cat(3,red,green,blue);
%         [M,N,K]=size(g);
%         I=reshape(g,M*N,3);
%         idx=find(epidermis_mask1);
%         I=double(I(idx,1:3));
%         [C,m]=covmatrix(I);
%         d=diag(C);
%         sd=sqrt(d);
%         for i=1:CC.NumObjects
%             mask1=false(size(initial_epidermis_mask));
%             mask1(CC.PixelIdxList{i})=1;
%             red=immultiply(mask1,IM(:,:,1));%??????
%             green=immultiply(mask1,IM(:,:,2));
%             blue=immultiply(mask1,IM(:,:,3));
%             f=cat(3,mean(red(:))*255,mean(green(:))*255,mean(blue(:))*255);
%             E30=colorseg('mahalanobis',f,round(min(sd))-10,m,C);
%             if E30==1
%                 epidermis_mask(CC.PixelIdxList{i})=1;
%             end
%         end
%     end
end
epidermis_mask=imclose(epidermis_mask,strel('disk',7));
epidermis_mask=imdilate(epidermis_mask,strel('disk',5));
epidermis_mask=imfill(epidermis_mask,'holes');
end