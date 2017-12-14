%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function is designed by Hongming Xu,
% Deptment of Eletrical and Computer Engineering,
% University of Alberta, Canada.  1th April, 2016
% If you have any problem feel free to contact me.
% Please address questions or comments to: mxu@ualberta.ca

% The techique is mainly based on the following paper:
% Xu, Hongming, et al. "Epidermis segmentation in skin histopathological images based on thickness measurement and k-means algorithm." (2015).

% Terms of use: You are free to copy,
% distribute, display, and use this work, under the following
% conditions. (1) You must give the original authors credit. (2) You may
% not use or redistribute this work for commercial purposes. (3) You may
% not alter, transform, or build upon this work. (4) For any reuse or
% distribution, you must make clear to others the license terms of this
% work. (5) Any of these conditions can be waived if you get permission
% from the authors.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function MainEpidermisSegmentation
close all;
clear all;

%% add the function into MATLAB searching path and enter the test dataset
p = mfilename('fullpath');
t = findstr(p,'\');
p = p(1:t(end));
addpath(p);
cd(p);
cd('WSI');
List = dir('*.tif');

%% coarse-to-fine epidermis segmentation
for i=4:length(List)
    imageName = List(i).name;
    IM = imread(imageName);
    
    maskEpidermis_PyraimdTop=XSegRoughEpidermis4ToBeEnhanced(IM,'R',0);                  % coarse segmentaton
    
    imagethick=XThicknessCal(maskEpidermis_PyraimdTop);                                  % thickness measurement
    
    [epidermis_mask,imagethick]=XSegEpidermis(IM,imagethick,maskEpidermis_PyraimdTop);   % fine segmentation
    
    figure,imshow(IM);
    B=bwboundaries(epidermis_mask);
    if ~isempty(B)
        boundary = B{1};
        hold on, plot(boundary(:,2), boundary(:,1), 'c', 'LineWidth', 2);
    end
    
end
