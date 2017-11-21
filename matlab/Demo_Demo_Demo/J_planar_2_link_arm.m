% Estimate end-effector position Jacobian of (horizontal planar) 3-DOF rotary joint arm.
%
%  J = J_planar_2_link_arm ( q, model )
%
%  in:
%      q       - joint angles
%      model   - model structure
%
%  out:
%      J       - end-effector position Jacobian
%
function J = J_planar_2_link_arm ( q, model )

L = model.L;
J(1,1,:) = -L(1)*sin(q(1,:))-L(2)*sin(q(1,:)+q(2,:));
J(1,2,:) =                  -L(2)*sin(q(1,:)+q(2,:));
J(2,1,:) =  L(1)*cos(q(1,:))+L(2)*cos(q(1,:)+q(2,:));
J(2,2,:) =                   L(2)*cos(q(1,:)+q(2,:));
J(3,1,:) = 1;
J(3,2,:) = 1;
