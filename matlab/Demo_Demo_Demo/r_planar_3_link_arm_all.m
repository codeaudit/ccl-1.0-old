function r = r_planar_3_link_arm_all ( Q, model )

L = model.L;
r1 = [L(1)*cos(Q(1));
    L(1)*sin(Q(1))];
r2 = [r1(1) + L(2)*cos(Q(1)+Q(2));
    r1(2) + L(2)*sin(Q(1)+Q(2))];
r3 = [r2(1)+L(3)*cos(Q(1)+Q(2)+Q(3));
      r2(2)+L(3)*sin(Q(1)+Q(2)+Q(3))];
r = [r1,r2,r3];