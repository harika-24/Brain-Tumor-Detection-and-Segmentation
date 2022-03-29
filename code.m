


### TO CONVERT THE PICTURE INTO GRAYSCALE:

RGB = imread('peppers.png');
imshow(RGB)
I = rgb2gray(RGB);
figure
imshow(I)


### MEDIAN FILTER 


I = imread('eight.tif');
figure, imshow(I)
K = medfilt2(I);


### Sobel METHOD
close all;
clc;
img=imread('a.jpg');
B=rgb2gray(img);
subplot(2,2,1)
imshow(B)
pause(2)
I=double(B);

for i=1:size(I,1)-2
for j=1:size(I,2)-2
%Sobel mask for x-direction:
mx=((2*I(i+2,j+1)+I(i+2,j)+I(i+2,j+2))-(2*I(i,j+1)+I(i,j)+I(i,j+2)));
%Sobel mask for y-direction:
my=((2*I(i+1,j+2)+I(i,j+2)+I(i+2,j+2))-(2*I(i+1,j)+I(i,j)+I(i+2,j)));

B(i,j)=sqrt(mx.^2+my.^2);
end
end
subplot(2,2,2)
imshow(B); title('Sobel gradient');
pause(2)
%Define a threshold value


Thresh=100;
B=max(B,Thresh);
B(B==round(Thresh))=0;
B=uint8(B);

subplot(2,2,3)
imshow(~B);title('Edge detected Image');




### SOBEL FILTER:

close all;
clc;
img=imread('a.jpg');
B=rgb2gray(img);
subplot(2,2,1)
imshow(B)
pause(2)
I=double(B);

for i=1:size(I,1)-2
for j=1:size(I,2)-2
%Sobel mask for x-direction:
mx=((2*I(i+2,j+1)+I(i+2,j)+I(i+2,j+2))-(2*I(i,j+1)+I(i,j)+I(i,j+2)));
%Sobel mask for y-direction:
my=((2*I(i+1,j+2)+I(i,j+2)+I(i+2,j+2))-(2*I(i+1,j)+I(i,j)+I(i+2,j)));

B(i,j)=sqrt(mx.^2+my.^2);
end
end
subplot(2,2,2)
imshow(B); title('Sobel gradient');
pause(2)
%Define a threshold value


Thresh=100;
B=max(B,Thresh);
B(B==round(Thresh))=0;
B=uint8(B);

subplot(2,2,3)
imshow(~B);title('Edge detected Image');




### REGION GROWING:

I = im2double(imread('Image_To_Read.tiff'));
figure, imshow(I)
%imtool(I);
Isizes = size(I); %size of the image

threshI = multithresh(I, 3); %thresholding for three regions

[m, n]=ginput(1); %pick one pixel of the region to be segmented
c = impixel(I, m, n); %value of the pixel picked
currPix = c(1); %current pixel

surr = [-1 0; 1 0; 0 -1; 0 1]; %create a mask which represents the four surrounding pixels
mem = zeros(Isizes(1)*Isizes(2), 3); %create a register to put the pixel coordinates and pixel value 
mem(1, :) = [m, n, currPix]; %insert initial picked pixel to the register
regSize = 1; %initial size

J = zeros(Isizes(1), Isizes(2)); %create another black image with the same size as the original image

init = 1;
posInList = 1;
k=1; %create the initial condition to run the loop
%The region growing algorithm.
while(k==1)
    
  for l=init:posInList %first pointer on the register
  for j=1:4 %second pointer for the neighboring pixels
        m1 = m + surr(j,1);
        n1 = n + surr(j,2);
        
        check=(m1>=1)&&(n1>=1)&&(m1<=Isizes(1))&&(n1<=Isizes(2)); %simple check if pixel position still inside the image
        
        current = impixel(I, m1, n1);
        currPix = current(1);
        if(check && currPix<=threshI(2) && (J(m1, n1)==0)) %check if it belongs to the thresholding boundary and if not set yet on the image we want to recreate
            posInList = posInList+1;
            mem(posInList, :) = [m1, n1, currPix]; %add the new pixel
            J(m1, n1) = 1;
        end
  end
  end
  if(posInList == init) %when there is no more pixels to add
      k = 0; %make k=0 to close the loop
  else
      init = init+1;
      m = mem(init, 1, :);
      n = mem(init, 2, :);
      k = 1; %keep running the loop
  end
end

imshow(J); %the segmented black and white region

