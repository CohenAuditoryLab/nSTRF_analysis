function [cleanedSTRF]=binarizeSTRFs(inputSTRF)
    [rows, cols] = size(inputSTRF);
    numPixels = double(rows * cols);
    threshold = 0.05 / numPixels;
    for i=1:rows
        for j=1:cols
            if(inputSTRF(i,j) > threshold)
                inputSTRF(i,j) = 1;
            else
                inputSTRF(i,j) = 0;
            end
        end
    end
    cleanedSTRF = inputSTRF;
end
