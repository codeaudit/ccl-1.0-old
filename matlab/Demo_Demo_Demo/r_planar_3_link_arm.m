% Estimate end-effector position of (horizontal planar) 3-DOF rotary joint arm.
%
%  r = r_planar_3_link_arm ( q, model )
%
%  in:
%      q       - joint angles
%      model   - model structure
%
%  out:
%      r       - end-effector position
%
function r = r_planar_3_link_arm ( q, model )

L = model.L;
r = r_planar_2_link_arm ( q, model ) + [ L(3)*cos(q(1,:)+q(2,:)+q(3,:));...
                                         L(3)*sin(q(1,:)+q(2,:)+q(3,:));...
                                         q(3,:)];

