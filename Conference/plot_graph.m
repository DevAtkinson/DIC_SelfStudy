% load('G_sin.mat');
load('G_final_sin.mat');
load('F_in_final.mat');
figure
[r,c]=size(F);
plot(1:1:c,F(20,:),'k')
hold on;
[r,c]=size(G_store{end});
plot(1:1:c,G_store{end}(20,:),'Color', [0.6 0.6 0.6],'Marker','*','MarkerEdgeColor','k','MarkerSize',4)

clear G_store
% load('G_normal.mat');
load('G_final_normal.mat');
[r,c]=size(G_store{end});
plot(1:1:c,G_store{end}(20,:),'Color', [0.6 0.6 0.6],'Marker','o','MarkerEdgeColor','k','MarkerSize',4)
plot(1:1:c,F(20,:),'k')

% clear G_store
% load('G_final_zero.mat');
% [r,c]=size(G_store{end})
% plot(1:1:c,G_store{end}(20,:),':')
legend('original subset','sin shape function','affine shape function')
title('Comparison of light intensity values over the subsets')
xlabel('Pixel [x-direction]')
ylabel('Light intensity')
axis([0 72 0 0.025])