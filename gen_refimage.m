clear;

original = double(rgb2gray(imread('cave01_00_flash.tif')))/255;
guide =  double(rgb2gray(imread('cave01_01_noflash.tif')))/255;
alpha = 15;
ws = 120;
eps = 0.01^2;
col_start = 120*5;

I = original(:,col_start+1:col_start+4*alpha+ws);
p = guide(:,col_start+1:col_start+4*alpha+ws);

%%
qred = guidedfilter(I, p, alpha, eps);
qred = (qred(:,2*alpha+1:end-2*alpha));
qred(qred>1) = 1;
qred(qred<0) = 0;
q = double(q)/255;
psnr(q, qred)
%%

q = guidedfilter(I, p, alpha, eps);

figure();
imshow([I, p, q], [0, 1]);
