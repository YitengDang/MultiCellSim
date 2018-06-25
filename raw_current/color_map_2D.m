rows = 500;
cols = 500;
image=zeros(rows, cols,3); %initialize
for row=1:rows
    for col=1:cols
        y = 1-(row-1)/(rows-1);
        b = 1-(col-1)/(cols-1);
        image(row, col, :) = [y y b];
    end
end
h=figure;
imshow(image)
set(gca, 'YDir', 'normal')
xlabel('Gene 1');
ylabel('Gene 2');
%set(gca, 'XTick', [0:10:500]);