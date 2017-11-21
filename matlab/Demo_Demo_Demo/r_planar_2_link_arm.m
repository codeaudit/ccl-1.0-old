% Estimate end-effector position of (horizontal planar) 3-DOF rotary joint arm.
%
%  r = r_planar_2_link_arm ( q, model )
%
%  in:
%      q       - joint angles
%      model   - model structure
%
%  out:
%      r       - end-effector position
%
function r = r_planar_2_link_arm ( q, model )

L = model.L;
r = [L(1)*cos(q(1,:)) + L(2)*cos(q(1,:)+q(2,:));...
     L(1)*sin(q(1,:)) + L(2)*sin(q(1,:)+q(2,:));...
     q(1)+q(2)];
 

