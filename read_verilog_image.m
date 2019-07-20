clear;
q = uint8(zeros( 1080, 1920));

for k = 1:16
    fprintf('Read Image %d\n', k)
    fileID = fopen(strcat('result/outputImage_stripe_',num2str(k-1)),'r');
    read_image = uint8(zeros(1080, 120));

    for i = 1:1080
        for j = 1:120
            read_image(i, j) = uint8(str2num(fgetl(fileID)));
        end
    end
   q(:,120*(k-1)+1:120*k) = read_image;
   %q = read_image;
    fclose(fileID);
    %typecast(uint8(hex2dec(f(1+8*(i-1):8+8*(i-1)))),'single');
end