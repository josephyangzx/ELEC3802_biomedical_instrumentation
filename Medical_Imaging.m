%% Lab 2: Medical imaging
% ELEC3802

clear;
clc;
clf;
%% I.Data
load mri;

% Set the parameters for plots
% Set the position 
x0=10;
y0=10;
% Set the height and width of the plot window
width=900;
height=500;

% axial (:,:,1); saggital (:,1,:); coronal (1,:,:)
% Axial view of slice 16
im_axial = squeeze(D(:,:,16));
figure(1)
set(gcf,'position',[x0,y0,width,height])
subplot(1,5,1:3);
imshow(im_axial, []);
title('Axial view at slice 16');
saveas(gcf,'axial_16','png');

% D (:,:,1) - Axial, squeeze(D(:,2,:))
% - Sagittal; squeeze (D(1,:,:))
% make E = squeeze(D); E(1,:) for coronal
% Sagittal view of slice 64
im_sagittal = squeeze(D(:,64,:));
subplot(1,5,4);
imshow(im_sagittal, []);
title('Sagittal view at slice 64');
saveas(gcf,'sagittal_64','png');

% Coronal view of slice 64
im_coronal = squeeze(D(64,:,:));
subplot(1,5,5);
imshow(im_coronal, []);
title('Coronal view at slice 64');
saveas(gcf,'coronal_64','png');
% Set the main title of the subplot
sgtitle('Anatomical View of Different Slices');
saveas(gcf,'sectionalViews','png');

%% II. Edge Filters

% Image gradient
% imggradientxy using prewitt kernel
[gradX,gradY] = imgradientxy(im_axial, 'prewitt');
% imggradient using prewitt kernel
[Gmag,Gdir] = imgradient(im_axial, 'prewitt');
% Convert the output to unit8
gradXU8 = uint8(gradX);
gradYU8 = uint8(gradY);
GmagU8 = uint8(Gmag);
GdirU8 = uint8(Gdir);
% Plot the outputs
figure();
subplot(2,2,1);
imshow(gradXU8, []);
title('Prewitt horizontal edges');
subplot(2,2,2);
imshow(gradYU8, []);
title('Prewitt vertical edges');
subplot(2,2,3);
imshow(GmagU8, []);
title('Prewitt Gmag');
subplot(2,2,4);
imshow(GdirU8, []);
title('Prewitt Gdir');
saveas(gcf,'prewitt','png');

% imggradientxy using sobel kernel
[gradXs,gradYs] = imgradientxy(im_axial, 'sobel');
% imggradient using sobel kernel
[GmagS,GdirS] = imgradient(im_axial, 'sobel');
% Convert the output to unit8
gradXsU8 = uint8(gradXs);
gradYsU8 = uint8(gradYs);
GmagSU8 = uint8(GmagS);
GdirSU8 = uint8(GdirS);
% Plot the oputputs
figure();
subplot(2,2,1);
imshow(gradXsU8, []);
title('Sobel horizontal edges');
subplot(2,2,2);
imshow(gradYsU8, []);
title('Sobel vertical edges');
subplot(2,2,3);
imshow(GmagSU8, []);
title('Sobel Gmag');
subplot(2,2,4);
imshow(GdirSU8, []);
title('Sobel Gdir');
saveas(gcf,'sobel','png');

% Canny edge default
[BW1,threshOut1] = edge(im_axial,'canny');
% Canny edge with prewitt
BW2 = edge(im_axial,'prewitt');
% Canny edge with sobel
BW3 = edge(im_axial,'sobel');
% Conver the output to unit8
BW1U8 = uint8(BW1);
BW2U8 = uint8(BW2);
BW3U8 = uint8(BW3);
% Plot the figures
figure();
subplot(2,1,1);
imshowpair(BW1U8,BW2U8,'montage');
title('Canny Prewitt');
subplot(2,1,2);
imshowpair(BW1U8,BW3U8,'montage');
title('Canny Sobel');
saveas(gcf,'canny_Edge','png');

% Change the threshold of Canny edge
[BW4,threshOut4] = edge(im_axial,'canny',0.3);
[BW5,threshOut5] = edge(im_axial,'canny',0.05);
% Conver the output to unit8
BW4U8 = uint8(BW4);
BW5U8 = uint8(BW5);
figure();
imshowpair(BW4U8,BW5U8,'montage');
title('Canny 0.3(L) & 0.05(R)');
saveas(gcf,'canny_Edge_threshold','png');

%% Kmeans clustering
% Ensure the imput is double
im_axial_double = double(im_axial);
% Get the size of im_axial_double
img__axial_size = size(im_axial_double);
% Proceed for k-means clustering
im_axial_reshape = reshape(im_axial_double, img__axial_size(1)*img__axial_size(2),1);

[idx1, C1] = kmeans(single(im_axial_reshape),4,'Display','iter');
imcluster1 = zeros(size(idx1));

for i = 1:4 
   imcluster1(find(idx1 == i))= C1(i); 
end

imcluster1 = reshape(imcluster1, img__axial_size(1),img__axial_size(2));

[idx2, C2] = kmeans(single(im_axial_reshape),8,'Display','iter');
imcluster2 = zeros(size(idx2));

for i = 1:8 
   imcluster2(find(idx2 == i))= C2(i); 
end

imcluster2 = reshape(imcluster2, img__axial_size(1),img__axial_size(2));

[idx3, C3] = kmeans(single(im_axial_reshape),20,'Display','iter');
imcluster3 = zeros(size(idx3));

for i = 1:20 
   imcluster3(find(idx3 == i))= C3(i); 
end

imcluster3 = reshape(imcluster3, img__axial_size(1),img__axial_size(2));

% Convert theoutput to unit8
imcluster1U8 = uint8(imcluster1);
imcluster2U8 = uint8(imcluster2);
imcluster3U8 = uint8(imcluster3);

% Set the dimension of the plots
x0=10;
y0=10;
width=2400;
height=470;
figure();
set(gcf,'position',[x0,y0,width,height])
subplot(1,3,1);
imagesc(imcluster1U8);
title('4 clusters (uint8)');
% saveas(gcf,'4_clusters_uint8','png');
subplot(1,3,2);
imagesc(imcluster2U8);
title('8 clusters (uint8)');
% saveas(gcf,'8_clusters_uint8','png');
subplot(1,3,3);
imagesc(imcluster3U8);
title('20 clusters (uint8)');
saveas(gcf,'clusters_uint8','png');

%% Plots in Double
figure();
imagesc(imcluster2);
title('8 clusters (double)');
saveas(gcf,'8_clusters','png');

figure();
imagesc(imcluster1);
title('4 clusters (double)');
saveas(gcf,'4_clusters','png');

figure();
imagesc(imcluster3);
title('20 clusters (double)');
saveas(gcf,'20_clusters','png');

%%
% Perform kmeans and show distance sum and number of iterations
[idx,C,sumD] = kmeans(double(im_axial_reshape),4, 'display', 'iter');

% The iter and sum are determine AFTER running the codes since there is no
% way to determine the number of interations
iter = 1:6;
sum = [485651 390208 325223 299100 289060 286506];
plot(iter, sum);
title('within-cluster sums of point-to-centroid distances VS # OF iteration');
saveas(gcf,'iter','png');

% END of Lab 2: Medical imaging%