function h = plot_arm_3link (h, Q, L, C)
% Plot_arm (Q, L, C)
% A function to plot arm
% Input:
%   Q               Joint space vector
%   L               Link lengths
%   C               Color
if ~exist('C', 'var')
    C = 'r' ;
end
r1 = zeros(2,1); % base

r2 = [L(1)*cos(Q(1));
    L(1)*sin(Q(1))];
r3 = [r2(1) + L(2)*cos(Q(1)+Q(2));
    r2(2) + L(2)*sin(Q(1)+Q(2))];
r4 = [r3(1)+L(3)*cos(Q(1)+Q(2)+Q(3));
      r3(2)+L(3)*sin(Q(1)+Q(2)+Q(3))];
h1 = plot([r1(1) r2(1)], [r1(2) r2(2)], 'LineStyle', '-', 'Color', C ) ;
h2 = plot([r2(1) r3(1)], [r2(2) r3(2)], 'LineStyle', '-', 'Color', C ) ;
h3 = plot([r3(1) r4(1)], [r3(2) r4(2)], 'LineStyle', '-', 'Color', C ) ;
h = [h1,h2,h3];
end
