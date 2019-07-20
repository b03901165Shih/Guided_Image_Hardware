clear;

% strip_ind: 0~15 for a full-HD image
stripe_ind = 2;


%addpath('He_code')
original = imread('cave01_00_flash.tif');
guide = imread('cave01_01_noflash.tif');
original = double(original(:,:,1))/255;
guide =  double(guide(:,:,1))/255;

alpha = 15;
ws = 120;
eps = 0.01^2;

original = padarray(original,[0 alpha*2],'both');
guide = padarray(guide,[0 alpha*2],'both');


col_start = 120*stripe_ind;

I = original(:,col_start+1:col_start+4*alpha+ws);
p = guide(:,col_start+1:col_start+4*alpha+ws);

[hei, wid] = size(I);
N = int64(boxfilter(ones(hei, wid), alpha));
N = N(:, 1+alpha:end-alpha);

mean_I = int64((boxfilter(I, alpha))*255) ;
mean_p = int64((boxfilter(p, alpha)) *255);
mean_Ip = int64((boxfilter(I.*p, alpha)) *255*255);
mean_II = int64((boxfilter(I.*I, alpha))*255*255);

% Remove first and last alpha columns
mean_I = mean_I(:, 1+alpha:end-alpha);
mean_p = mean_p(:, 1+alpha:end-alpha);
mean_Ip = mean_Ip(:, 1+alpha:end-alpha);
mean_II = mean_II(:, 1+alpha:end-alpha);

cov_Ip = mean_Ip.*N - mean_I .* mean_p; % this is the covariance of (I, p) in each local patch.
var_I = mean_II.*N - mean_I .* mean_I;

%cov_Ip = double(cov_Ip);
%var_I = double(var_I);

frac_bits = 8;
frac_bit_left_in_a = 6; 
frac_bit_left_in_b = 1;

a = int64(floor(double(cov_Ip.*2^frac_bits) ./ double(var_I + eps*255*255*(2*alpha+1)^4))); % Eqn. (5) in the paper;
mask = cov_Ip<0;
a = a+int64(mask); %add 1 if negative 

b_a = int64(floor(double(a)/double(2^(frac_bits-frac_bit_left_in_a))));
numerator = (mean_p.*2^frac_bit_left_in_a - b_a .* mean_I);
b = int64(floor(double(numerator)./double(N.*2^(frac_bit_left_in_a-frac_bit_left_in_b)))); % Only one fractional bit
mask = numerator<0;
b = b+int64(mask); %add 1 if negative 

a = double(a);
b = double(b);

mean_a = int64(boxfilter(a, alpha)) ;
mean_b = int64(boxfilter(b, alpha)) ;

% Remove first and last alpha columns
mean_a = mean_a(:, 1+alpha:end-alpha);
mean_b = mean_b(:, 1+alpha:end-alpha);
N = N(:, 1+alpha:end-alpha);

I = int64(I*255);
p = int64(p*255);

q = (mean_a .* I(:,1+2*alpha:end-2*alpha) + mean_b.*2^(frac_bits-frac_bit_left_in_b))./(N.*2^frac_bits); % Eqn. (8) in the paper;

q(q>255) = 255;
q(q<0) = 0;
q = uint8(q);

fprintf('Start writing pattern...\n');

I_path = 'I.dat';
P_path = 'p.dat';
ak_path = 'ak.dat';
bk_path = 'bk.dat';
q_path = 'q.dat';
write2File = true;

if(write2File)
    % Write file
    fileIDI = fopen(I_path,'w');
    fileIDP = fopen(P_path,'w');
    for i = 1:(size(I,1))
        for j = 1:size(I,2)
            fprintf(fileIDI,'%s\n',dec2bin(I(i,j), 8));
            fprintf(fileIDP,'%s\n',dec2bin(p(i,j), 8));
        end
    end
    fclose(fileIDI);
    fclose(fileIDP);

    fileID = fopen(q_path,'w');
    for i = 1:size(q,1)
        for j = 1:size(q,2)
            s = dec2bin(q(i,j), 8);
            fprintf(fileID,'%s\n',s);
        end
    end
    fclose(fileID);
    
    %Third Stage
    %{
    fileID = fopen(ak_path,'w');
    for i = 1:size(mean_a,1)
        for j = 1:size(mean_a,2)
            s = dec2bin(typecast(int32(mean_a(i,j)),'uint32'),32);
            fprintf(fileID,'%s\n',s(32-23+1:end));
        end
    end
    fclose(fileID);
    % Write file
    fileID = fopen(bk_path,'w');
    for i = 1:size(mean_b,1)
        for j = 1:size(mean_b,2)
            s = dec2bin(typecast(int32(mean_b(i,j)),'uint32'),32);
            fprintf(fileID,'%s\n',s(32-24+1:end));
        end
    end
    fclose(fileID);
    %}
end
    
    %Second Stage
    %{
    % Write file
    fileID = fopen(ak_path,'w');
    for i = 1:size(a,1)
        for j = 1:size(a,2)
            s = dec2bin(typecast(int16(a(i,j)),'uint16'),16);
            fprintf(fileID,'%s\n',s(16-13+1:end));
        end
    end
    fclose(fileID);
    % Write file
    fileID = fopen(bk_path,'w');
    for i = 1:size(a,1)
        for j = 1:size(a,2)
            s = dec2bin(typecast(int16(b(i,j)),'uint16'),16);
            fprintf(fileID,'%s\n',s(16-14+1:end));
        end
    end
    fclose(fileID);
end
%}


%{
if(write2File)
    % Write file
    fileID = fopen(input_path,'w');
    for i = 1:size(stripe,1)
        for j = 1:size(stripe,2)
            fprintf(fileID,'%s\n',dec2bin(stripe(i,j), 8));
        end
    end
    fclose(fileID);
    % Write file
    fileID = fopen(gold_path,'w');
    for i = 1:size(stripe,1)
        for j = 1:size(stripe,2)
            if(i==1)
                fprintf(fileID,'%s',dec2bin(0, 8));
            else
                fprintf(fileID,'%s',dec2bin(stripe(i-1,j), 8));
            end
            if(j<=2*alpha+1)
                fprintf(fileID,'%s\n',dec2bin(0, 8));
            else
                fprintf(fileID,'%s\n',dec2bin(stripe(i,j-2*alpha-1), 8));
            end
        end
    end
    fclose(fileID);
end
%}