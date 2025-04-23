clc, clear


tic
initial_conditions = [-0.03         % d_f1
                       0.00         % d_pg1
                      -0.00         % d_m1
                       0.00         % d_ptie1
                      -0.03         % d_f2
                       0.00         % d_pg2
                      -0.00];       % dpm2
% generating mesh for f1 and f2
x = -0.4:0.013:0.4; %f1
y = x; %f2
[X,Y] = meshgrid(x,y);

evaluation = fillGrid(X,Y, initial_conditions); 


% Visualize the results using a heatmap
figure;
imagesc(x, y, evaluation);
colormap([1 0.5 0.5; 0.8 1 0.5]); % 
xlabel('\Delta f_1');
ylabel('\Delta f_2');
fontname(gcf,"Garamond")
fontsize(gcf,20,"pixels")
set(gcf,'Position',[10 10 750 500])
title('Mesh grid, f1 and f2');
toc


function [res] = fillGrid(X, Y, initial_conditions)
    [rows,cols] = size(X);    
    res = zeros(rows, cols);

    for i = 1:rows
        for j = 1:cols
            try
                initial_conditions(1) = X(i,j);
                initial_conditions(5) = Y(i,j);
                RegulatorControllerMPC(initial_conditions)  
                res(i,j) = 1;
            catch exception
                res(i,j) = 0;
            end
        end
    end
end


