function visualise_move_dual_3link (figh, R_l,R_r, X_l,X_r, L,varargin)
% visualise_move (R, X, L)
% Visualisation of the arm movement
% Input:
%   R               Task space movements
%   X               Joint space movements
%   L               Arm link length
%   varargin        Arguments for figure title and position
dim_n = length(R_l) ;
xmin = -1.5 ; xmax = 2 ;
ymin = -2.5 ; ymax = 2.5 ;

fig_handle = figure(figh); hold on;
if nargin > 3
    title_txt = varargin{1};title(title_txt);
    pos       = varargin{2};set(fig_handle, 'Position', pos);
end
% axis equal
xlim([xmin,xmax]) ; ylim([ymin,ymax]) ;
R1_l = squeeze([R_l(1,1,:); R_l(2,1,:)]);
R2_l = squeeze([R_l(1,2,:); R_l(2,2,:)]);
R3_l = squeeze([R_l(1,3,:); R_l(2,3,:)]);

plot( R1_l(1,:),R1_l(2,:), '--','LineWidth', 2) ;
plot( R2_l(1,:),R2_l(2,:),'--', 'LineWidth', 2) ;
plot( R3_l(1,:),R3_l(2,:),'--', 'LineWidth', 2) ;

R1_r = squeeze([R_r(1,1,:); R_r(2,1,:)]);
R2_r = squeeze([R_r(1,2,:); R_r(2,2,:)]);
R3_r = squeeze([R_r(1,3,:); R_r(2,3,:)]);

plot( R1_r(1,:),R1_r(2,:), '--','LineWidth', 2) ;
plot( R2_r(1,:),R2_r(2,:),'--', 'LineWidth', 2) ;
plot( R3_r(1,:),R3_r(2,:),'--', 'LineWidth', 2) ;


xlabel('x');
ylabel('y');
% stroboscopic plot of arm
cc = rand;
c = [cc,cc,cc] ;
h1 = [];
h1 = plot_arm_3link (h1,X_l(:,1), L, c) ;axis tight;pause (0.1) ;
h2= [];
h2 = plot_arm_3link (h2,X_r(:,1), L, c) ;axis tight;pause (0.1) ;
c = [cc,0,0] ;
h3 = [];h4 = [];
for i=2:dim_n
    h3 = plot_arm_3link (h3,X_l(:,i), L, c) ;
    h4 = plot_arm_3link (h4,X_r(:,i), L, c) ;
    if i ~=dim_n
        pause (0.1) ;
        h3.reset; 
        h4.reset; 
    end
%     axis tight;
end
end
