%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% this is a function used to calculate the thickness of the epidermis.

% Input:
%   -IM    RGB image
%   -maskimage  the epidermis region
% Output:
%   -maskEpidermis    a logical matrix indicate the position of the
%                   epidermis
% Key Threshold:
%   -TAxisRatio % we define the enlonged one as the AxisRatio > TAxisRatio
%   -TAreaofROI the threshold for the size of area that we think is noise

% (c) Edited by Hongming Xu,
% Deptment of Eletrical and Computer Engineering,
% University of Alberta, Canada.  20th Oct, 2013
% If you have any problem feel free to contact me.
% Please address questions or comments to: mxu@ualberta.ca
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function imagethick=XThicknessCal(maskimage,shown)
if ~exist('shown','var')
    shown=0;
end

% maskimage=imfill(maskimage,'holes');
%% The first step -- smoothen the binary image
se=strel('disk',10);
maskimage_open=imclose(maskimage,se);
% maskimage_open=imfill(maskimage_open,'holes');

% h=fspecial('gaussian',17,10);
% maskimage_open=imfilter(maskimage,h,'replicate');
maskimage_open=imfill(maskimage_open,'holes');
CC=bwconncomp(maskimage_open);
numPixels=cellfun(@numel,CC.PixelIdxList);   %% find the maximum component
[biggest,idx]=max(numPixels);
maskimage_smooth=zeros(size(maskimage));
maskimage_smooth(CC.PixelIdxList{idx})=1;

%% The second step -- calculate the skeleton of the epidermis region
s2=bwmorph(maskimage_smooth,'thin',Inf);
% s2=bwmorph(maskimage_smooth,'skel',Inf);
if shown
    [r1,c1]=find(s2==1);
    show(maskimage_smooth,12);
    hold on,plot(c1,r1,'r.');
end     

%% The third step -- calculate the usable end points
s3=bwmorph(s2,'endpoints');

[r3,c3]=find(s3==1);
X=[r3,c3];Y=[r3,c3];
D=sqrt(sum(abs(repmat(permute(X,[1 3 2]),[1 length(Y) 1])...
    -repmat(permute(Y,[3 1 2]),[length(X) 1 1])).^2,3));
sortedValues=unique(D(:));
if length(sortedValues)>10            % select the top 10 pairs
    bigValues=sortedValues(end-9:end);% 10 big values
    bigIndex=ismember(D,bigValues);
else
    bigIndex=ismember(D,sortedValues);
end
Lindex=tril(bigIndex);
[rendpoints,cendpoints]=find(Lindex==1);

if shown
    [r3,c3]=find(s3==1);
    hold on,plot(c3,r3,'g+','MarkerSize',14,'LineWidth',3);
    
    % selected end points
    eindex=unique([rendpoints;cendpoints]);
    hold on,plot(c3(eindex),r3(eindex),'g+','MarkerSize',14,'LineWidth',3);hold off;
end

%% The fourth step -- calculate the skeleton of region
lmax=0;
for i=1:length(rendpoints)
    ip1=rendpoints(i);ip2=cendpoints(i);
    D1=bwdistgeodesic(s2,c3(ip1),r3(ip1),'quasi-euclidean');
    D2=bwdistgeodesic(s2,c3(ip2),r3(ip2),'quasi-euclidean');
    D=D1+D2;
    D=round(D*8)/8;
    D(isnan(D))=inf;
    minvalue=min(D(:));
    th1=minvalue+0.5;   % error with 0.5 is accepted
    [rt,ct]=find(D<th1);
    
    if length(rt)>lmax
        rf=rt;cf=ct;lmax=length(rt);
        
        t1=ip1;t2=ip2; %% add for drawing figure
    end
    if shown
        skeleton_path=zeros(size(maskimage));
        linearInd=sub2ind(size(maskimage),rt,ct);
        skeleton_path(linearInd)=1;
        P = imoverlay(s2, imdilate(skeleton_path, ones(3,3)), [1 0 0]);
        show(P,4);
        hold on,plot(c3(ip1),r3(ip1),'g.','MarkerSize',15);
        plot(c3(ip2),r3(ip2),'g.','MarkerSize',15);
        hold off;
    end
end

if shown
    skeleton_path=zeros(size(maskimage));
    linearInd=sub2ind(size(maskimage),rf,cf);
    skeleton_path(linearInd)=1;
    P = imoverlay(maskimage_smooth, imdilate(skeleton_path, ones(8,8)), [1 0 0]);
    show(P,5);
     hold on,plot(c3(t1),r3(t1),'g+','MarkerSize',14,'LineWidth',3);
        plot(c3(t2),r3(t2),'g+','MarkerSize',14,'LineWidth',3);
end



% lmax=0; 
% for i=1:length(r3)-1
%     for j=i+1:length(r3)
%         D1=bwdistgeodesic(s2,c3(i),r3(i),'quasi-euclidean');
%         D2=bwdistgeodesic(s2,c3(j),r3(j),'quasi-euclidean');
%         D=D1+D2;
%         D=round(D*8)/8;
%         D(isnan(D))=inf;
%         minvalue=min(D(:));
%         th1=minvalue+0.5;   % error with 0.5 is accepted
%         [rt,ct]=find(D<th1);
%         %              skeleton_path=imregionalmin(D);
%         if length(rt)>lmax
%             rf=rt;cf=ct;lmax=length(rt);
%         end
% %         if shown
% %             skeleton_path=zeros(size(maskimage));
% %             linearInd=sub2ind(size(maskimage),rt,ct);
% %             skeleton_path(linearInd)=1;
% %             P = imoverlay(s2, imdilate(skeleton_path, ones(3,3)), [1 0 0]);
% %             show(P,4);
% %             hold on,plot(c3(i),r3(i),'g.','MarkerSize',15);
% %             plot(c3(j),r3(j),'g.','MarkerSize',15);
% %             hold off;
% %         end
%     end
% end



%% The fifth step -- calculate the boundaries points and down sample the skeleton
curB=[];
B=bwboundaries(maskimage_smooth,'noholes');
for i=1:length(B)
    cur=B{i};
    curB=[curB;cur];
end

% skeleton_points=[rf,cf];
% dsample_index=downsample(1:length(skeleton_points),10);
% dsample_points=skeleton_points(dsample_index,:);

temp2=20; % down-sampling factor
skeleton_path=zeros(size(maskimage));
linearInd=sub2ind(size(maskimage),rf,cf);  
skeleton_path(linearInd)=1;
sp1=bwmorph(skeleton_path,'endpoints');
[r5,c5]=find(sp1==1);
boundaries=bwtraceboundary(skeleton_path,[r5(1),c5(1)],'N');
skeleton_points=boundaries(1:round(length(boundaries)/2),:);

skeleton_points_filter=XFilter(skeleton_points,201);
dsample_points=downsample(skeleton_points_filter,temp2);

if shown
    show(maskimage_smooth,13);
    hold on,plot(skeleton_points(:,2),skeleton_points(:,1),'r.');
    plot(dsample_points(:,2),dsample_points(:,1),'b*'); hold off;
end

%% The sisth step --- calculate the thickness
if shown
    bw=zeros(size(maskimage_smooth));
    imshow(bw,13);
    hold on,plot(curB(:,2),curB(:,1),'g.');
%     hold on,plot(skeleton_points_filter(:,2),skeleton_points_filter(:,1),'r.');
    hold on,plot(dsample_points(:,2),dsample_points(:,1),'y*');
end
Dis=[];

for i=2:length(dsample_points)-1
    % a-- calculate the slope the of the vertical line
    
     p1=skeleton_points_filter((i-1)*temp2,:);
     p2=skeleton_points_filter((i-1)*temp2+2,:);
%     p11=skeleton_points((i-1)*20-1,:);
%     p12=skeleton_points((i-1)*20,:);
%     p21=skeleton_points((i-1)*20+2,:);
%     p22=skeleton_points((i-1)*20+3,:);
%     p1=(p11+p12)/2;p2=(p21+p22)/2;
    if p1(1,2)==p2(1,2)
        slope=0; % vertical line
    else
        if p1(1,1)==p2(1,1)
            slope=Inf;
        else
            k1=(p2(1,1)-p1(1,1))/(p2(1,2)-p1(1,2));
            k2=-1/k1;
            slope=k2; 
%             BB=dsample_points(i,1)-k2*dsample_points(i,2);
        end
    end
    % b-- find the intersecton points
    tempd=repmat(dsample_points(i,:),length(curB),1);
    dis_temp=(curB(:,1)-tempd(:,1))./(curB(:,2)-tempd(:,2));
%     if slope==Inf;
%         ind3=find(dis_temp==Inf);
%     else
        dif_temp=atan(dis_temp)-atan(slope);
        ind3=find((abs(dif_temp)<0.05)); % error within 0.05 is accepted
%     end
    
    %% determine the side of a point lies
    A=p2(1,1)-p1(1,1);B=p1(1,2)-p2(1,2);C=p2(1,2)*p1(1,1)-p1(1,2)*p2(1,1);
    curp=curB(ind3,:);
    D=curp(:,2).*A+curp(:,1).*B+C;
    indr=find(D<0);
    indl=find(D>0);
    if length(indr)>0 & length(indl)>0
        d1=sqrt(sum(abs(curp(indr,:)-repmat(dsample_points(i,:),size(curp(indr,:),1),1)).^2,2));
        [mvalue1,ivalue1]=min(d1);
        
        d2=sqrt(sum(abs(curp(indl,:)-repmat(dsample_points(i,:),size(curp(indl,:),1),1)).^2,2));
        [mvalue2,ivalue2]=min(d2);
        
        
        r11=round((dsample_points(i,1)+curp(indr(ivalue1),1))/2);
        c11=round((dsample_points(i,2)+curp(indr(ivalue1),2))/2);
        
        
        r22=round((dsample_points(i,1)+curp(indl(ivalue2),1))/2);
        c22=round((dsample_points(i,2)+curp(indl(ivalue2),2))/2);
        
        IN=inpolygon([r11 r22],[c11 c22],curB(:,1),curB(:,2));
        
        if sum(IN)==2
            Dis=[Dis,mvalue1+mvalue2];
%                   hold on,plot([dsample_points(i,2),curp(indr(ivalue1),2)],[dsample_points(i,1),curp(indr(ivalue1),1)],'r','LineWidth',3);
%                   hold on,plot([dsample_points(i,2),curp(indl(ivalue2),2)],[dsample_points(i,1),curp(indl(ivalue2),1)],'w','LineWidth',3);
%                   hold on,plot([,curp(indl(ivalue2),2),curp(indr(ivalue1),2)],[curp(indl(ivalue2),1),curp(indr(ivalue1),1)],'r','LineWidth',3);
        end
    end
    
end
imagethick=mean(Dis);
% Vthick=var(Dis);
end