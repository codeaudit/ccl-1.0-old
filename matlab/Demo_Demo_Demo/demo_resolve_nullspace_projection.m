function demo_resolve_nullspace_projection
rng('default');
settings.dim_x = [6];
settings.dim_u = [6];
settings.dim_r = [4];
settings.dt = [0.01];
settings.N = [20];
[X,Y,Pi] = data_generation(settings);
model.L = [1,1,1];
%% For left
options.dim_b = 20;
options.dim_r = 2;
phi_= [1 0 0;
       0 1 0];
J_l  = @(q)J_planar_3_link_arm(q,model);
PHI_l = @(q)(phi_*J_l(q));
[optimal_l] = learn_lambda_ccl (Y(1:3,:), X(1:3,:), PHI_l, options);
%% read C code model
optimal_cl = optimal_l;
optimal_cl.w = [];
[optimal_cl.c,optimal_cl.s2,optimal_cl.w] = read_to_matlab('learn_lambda_model_l.txt',optimal_cl.dim_x,optimal_cl.dim_k,optimal_cl.dim_r,optimal_cl.dim_b);
optimal_cl.phi       = @(x) phi_gaussian_rbf ( x, optimal_cl.c,  optimal_cl.s2 ); 
optimal_cl.lambda  = @(q) predict_lambda (q, optimal_cl, eye(optimal_cl.dim_r)) ;

%% For right
options.dim_b = 20;
options.dim_r = 2;
phi_= [1 0 0;
       0 1 0];
J_r  = @(q)J_planar_3_link_arm(q,model);
PHI_r = @(q)(phi_*J_r(q));
[optimal_r] = learn_lambda_ccl (Y(4:6,:), X(4:6,:), PHI_r, options);
%% read C code model
optimal_cr = optimal_r;
optimal_cr.w = [];
[optimal_cr.c,optimal_cr.s2,optimal_cr.w] = read_to_matlab('learn_lambda_model_r.txt',optimal_cr.dim_x,optimal_cr.dim_k,optimal_cr.dim_r,optimal_cr.dim_b);
optimal_cr.phi       = @(x) phi_gaussian_rbf ( x, optimal_cr.c,  optimal_cr.s2 ); 
optimal_cr.lambda  = @(q) predict_lambda (q, optimal_cr, eye(optimal_cr.dim_r)) ;

optimal.lambda = @(q)([optimal_l.lambda(q(1:3)) 0 0;0 0 optimal_r.lambda(q(4:6))]);
%% read C code model
optimal_c.lambda = @(q)([optimal_cl.lambda(q(1:3)) 0 0;0 0 optimal_cr.lambda(q(4:6))]);
% %% evaluate
[nPPE, vPPE, uPPE] = get_ppe_alpha (optimal_r.f_proj, X(4:6,:), Pi(4:6,:), Y(4:6,:));
%% Ground true nullspace projection
lambda = @(q)env(q);
A = @(q)(lambda(q)*PHI(q));
%% reproduction
settings.N = 2;
for i = 1:settings.N
    q0_ = [-2.09 0 0 -1.04 0 0]';
    rd  = [-0.3*rand-1,-2,0.3*rand+1,-2]';
    q0 = go_start_pose(q0_,rd);
    rn = rand;
    title = 'Reproductions uisng True constraint';
    data_reproduction(settings,q0,rn,lambda,title); hold on; pause;
    title = 'Reproductions uisng Learned constraint';
    data_reproduction(settings,q0,rn,optimal.lambda,title);pause;
    title = 'Reproductions uisng Learned constaint in C code';
    data_reproduction(settings,q0,rn,optimal_c.lambda,title);hold on;%% C code model
    pause;
    clf;
end
close all;
end

function  [X,Y,Pi] = data_generation(settings)
N = settings.N;
dt = settings.dt;
model.L = [1,1,1];
fk_l = @(q)r_planar_3_link_arm(q(1:3),model);
fk_all_l = @(q)r_planar_3_link_arm_all(q(1:3),model);
J_l  = @(q)J_planar_3_link_arm(q(1:3),model);

fk_r = @(q)r_planar_3_link_arm(q(4:6),model);
fk_all_r = @(q)r_planar_3_link_arm_all(q(4:6),model);
J_r  = @(q)J_planar_3_link_arm(q(4:6),model);
lambda = @(q)env(q);
J = @(q)([J_l(q), zeros(3,3);zeros(3,3),J_r(q)]);
phi_= [1 0 0 0 0 0;
    0 1 0 0 0 0;
    0 0 0 1 0 0
    0 0 0 0 1 0];
Phi = @(q)(phi_*J(q));

X = [];Y = [];Pi = [];
q0_ = [-2.09 0 0 -1.04 0 0]';
for i = 1:N
    rd  = [-0.3*rand-1,-2,0.3*rand+1,-2]';
    q0 = go_start_pose(q0_,rd);
    X_i = []; Y_i = []; R_li = [];R_ri = [];Pi_ = [];
    n = 1;
    q = q0;
    r = zeros(6,50);
    beta = 0.3*rand+0.5;
    while n < 51
        r = [fk_l(q);fk_r(q)];
        q_null_target = [];
        if r(1) <= -0.5
            q_null_target = [-pi/2,rand+1.2,rand+1.2];
        elseif r(1) > -0.5
            q_null_target = [-pi,rand+1.2,rand+1.2];
        end
        if r(4) >=0.5
            q_null_target = [q_null_target, -pi/2,-rand-1.2,-rand-1.2];
        elseif r(4)< 0.5
            q_null_target = [q_null_target, 0,-rand-1.2,-rand-1.2];
        end   
%         if r(1)<= -0.5 && r(4)>= 0.5
%             q_null_target = [-pi/2,rand+1.2,rand+1.2,-pi/2,-rand-1.2,-rand-1.2]';
%         elseif  r(1)>= -0.5 && r(4)<= 0.5
%             q_null_target = [-pi,rand+1.2,rand+1.2,0,-rand-1.2,-rand-1.2]';
%         elseif  r(1)> -0.5 && r(4)> 0.5
%             q_null_target = [-pi,rand+1.2,rand+1.2,-pi/2,-rand-1.2,-rand-1.2]';
%         elseif  r(1)< -0.5 && r(4)< 0.5
%             q_null_target = [-pi/2,rand+1.2,rand+1.2,0,-rand-1.2,-rand-1.2]';
%         end
        
        pi_ = @(q)(beta*(q_null_target'-q));
        A = lambda(q)*Phi(q);
        pinvA = pinv(A);
        N = eye(6) - pinvA*A;
        u_ts = zeros(6,1);
        u_ns = N*pi_(q);
        u =  u_ts+ u_ns;
        X_i = [X_i q];
        Y_i = [Y_i u];
        Pi_ = [Pi_,pi_(q)];
        R_li(:,:,n) = fk_all_l(q);
        R_ri(:,:,n) = fk_all_r(q);
        q = q + dt*u;
        n = n + 1;
    end
%     if i <= 3
%         fig = figure(1);
%         visualise_move_dual_3link(fig,R_li,R_ri,X_i(1:3,:),X_i(4:6,:),model.L,'Resolving nullspace projection -- a Grasping task',[0,500,600,600]);
%     end
    X = [X,X_i];
    Y = [Y,Y_i];
    Pi = [Pi,Pi_];
end
end

function  data_reproduction(settings,q0,rn,lambda,tt)
dt = settings.dt;
model.L = [1,1,1];
fk_l = @(q)r_planar_3_link_arm(q(1:3),model);
fk_all_l = @(q)r_planar_3_link_arm_all(q(1:3),model);
J_l  = @(q)J_planar_3_link_arm(q(1:3),model);

fk_r = @(q)r_planar_3_link_arm(q(4:6),model);
fk_all_r = @(q)r_planar_3_link_arm_all(q(4:6),model);
J_r  = @(q)J_planar_3_link_arm(q(4:6),model);
J = @(q)([J_l(q), zeros(3,3);zeros(3,3),J_r(q)]);
phi_= [1 0 0 0 0 0;
    0 1 0 0 0 0;
    0 0 0 1 0 0
    0 0 0 0 1 0];
Phi = @(q)(phi_*J(q));
X_i = []; Y_i = []; R_li = [];R_ri = [];
n = 1;
q = q0;
beta = 0.5;
while n < 50
    r = [fk_l(q);fk_r(q)];
    q_null_target = [];
    if r(1) <= -0.5
        q_null_target = [-pi/2,rn+1.2,rn+1.2];
    elseif r(1) > -0.5
        q_null_target = [-pi,rn+1.2,rn+1.2];
    end
    if r(4) >=0.5
        q_null_target = [q_null_target, -pi/2,-rn-1.2,-rn-1.2];
    elseif r(4)< 0.5
        q_null_target = [q_null_target, 0,-rn-1.2,-rn-1.2];
    end
%     if r(1)<= -0.5 && r(4)>= 0.5
%         q_null_target = [-pi/2,rn+1.2,rn+1.2,-pi/2,-rn-1.2,-rn-1.2]';
%     elseif  r(1)>= -0.5 && r(4)<= 0.5
%         q_null_target = [-pi,rn+1.2,rn+1.2,0,-rn-1.2,-rn-1.2]';
%     elseif  r(1)> -0.5 && r(4)> 0.5
%         q_null_target = [-pi,rn+1.2,rn+1.2,-pi/2,-rn-1.2,-rn-1.2]';
%     elseif  r(1)< -0.5 && r(4)< 0.5
%         q_null_target = [-pi/2,rn+1.2,rn+1.2,0,-rn-1.2,-rn-1.2]';
%     end
    pi_ = @(q)(beta*(q_null_target'-q));
    A = lambda(q)*Phi(q);
    pinvA = pinv(A);
    N = eye(6) - pinvA*A;
    u_ts = zeros(6,1);
    u_ns = N*pi_(q);
    u =  u_ts+ u_ns;
    X_i = [X_i q];
    Y_i = [Y_i u];
    R_li(:,:,n) = fk_all_l(q);
    R_ri(:,:,n) = fk_all_r(q);
    q = q + dt*u;
    n = n + 1;
end
fig = figure(5);
visualise_move_dual_3link(fig,R_li,R_ri,X_i(1:3,:),X_i(4:6,:),model.L,tt,[0,500,600,600]);
end

function  data_reproduction_alpha(settings,q0,rn,A)
dt = settings.dt;
model.L = [1,1,1];
fk_l = @(q)r_planar_3_link_arm(q(1:3),model);
fk_all_l = @(q)r_planar_3_link_arm_all(q(1:3),model);
J_l  = @(q)J_planar_3_link_arm(q(1:3),model);

fk_r = @(q)r_planar_3_link_arm(q(4:6),model);
fk_all_r = @(q)r_planar_3_link_arm_all(q(4:6),model);
J_r  = @(q)J_planar_3_link_arm(q(4:6),model);
J = @(q)([J_l(q), zeros(3,3);zeros(3,3),J_r(q)]);
phi_= [1 0 0 0 0 0;
    0 1 0 0 0 0;
    0 0 0 1 0 0
    0 0 0 0 1 0];
Phi = @(q)(phi_*J(q));
X_i = []; Y_i = []; R_li = [];R_ri = [];
n = 1;
q = q0;
r = [fk_l(q);fk_r(q)];
while n < 50
    if r(1)<= -0.5 && r(4)>= 0.5
        q_null_target = [-pi/2,rn+1.2,rn+1.2,-pi/2,-rn-1.2,-rn-1.2]';
    elseif  r(1)>= -0.5 && r(4)<= 0.5
        q_null_target = [-pi,rn+1.2,rn+1.2,0,-rn-1.2,-rn-1.2]';
    elseif  r(1)> -0.5 && r(4)> 0.5
        q_null_target = [-pi,rn+1.2,rn+1.2,-pi/2,-rn-1.2,-rn-1.2]';
    elseif  r(1)< -0.5 && r(4)< 0.5
        q_null_target = [-pi/2,rn+1.2,rn+1.2,0,-rn-1.2,-rn-1.2]';
    end
    pi_ = @(q)(1*(q_null_target-q));
    pinvA = pinv(A(q));
    N = eye(6) - pinvA*A(q);
    u_ts = zeros(6,1);
    u_ns = N*pi_(q);
    u =  u_ts+ u_ns;
    X_i = [X_i q];
    Y_i = [Y_i u];
    R_li(:,:,n) = fk_all_l(q);
    R_ri(:,:,n) = fk_all_r(q);
    q = q + dt*u;
    r = [fk_l(q);fk_r(q)];
    n = n + 1;
end
fig = figure(5);
visualise_move_dual_3link(fig,R_li,R_ri,X_i(1:3,:),X_i(4:6,:),model.L,'Resolving nullspace projection -- a Grasping task',[0,500,600,600]);
end

function lambda = env(x)
model.L = [1,1,1];
fk_l = @(x)r_planar_3_link_arm(x(1:3),model);
fk_r = @(x)r_planar_3_link_arm(x(4:6),model);
r = [fk_l(x);fk_r(x)];
if r(1)<=-0.5
    lambda = [0 1 0 0];
elseif r(1)> -0.5
    lambda = [1 0 0 0;
              0 1 0 0];
end
if r(4)>=0.5
    lambda = [lambda;
              0 0 0 1];
elseif r(4) < 0.5
    lambda = [lambda;
              0 0 1 0;
              0 0 0 1];
end
% if r(1)<= -0.5 && r(4)>=0.5
%     lambda = [0 1 0 0;
%         0 0 0 1];
% elseif r(1)>= -0.5 && r(4)<=0.5
%     lambda = eye(4);%lambda = [1 0 0 0;
%        % 0 0 1 0];
% elseif r(1)< -0.5 && r(4)<0.5
%     lambda = eye(4);%lambda = [0 1 0 0;
%         %0 0 1 0];
% elseif r(1)> -0.5 && r(4)>0.5
%     lambda = eye(4);%lambda = [1 0 0 0;
%         %0 0 0 1];
% end
end

function q0 = go_start_pose(q0_,rd)
dt = 0.01;
model.L = [1,1,1];
lambda = [1 0 0 0 0 0;
    0 1 0 0 0 0;
    0 0 0 1 0 0;
    0 0 0 0 1 0];
q_default = [-2.09 0.52 0.52 -1.04 -0.52 -0.52]'+0.1*rand;
fk_l = @(q)r_planar_3_link_arm(q(1:3),model);
fk_r = @(q)r_planar_3_link_arm(q(4:6),model);
J_l  = @(q)J_planar_3_link_arm(q(1:3),model);
J_r  = @(q)J_planar_3_link_arm(q(4:6),model);
J = @(q)(lambda*[J_l(q) zeros(3,3);zeros(3,3),J_r(q)]);
r = @(q)([fk_l(q);fk_r(q)]);
q = q0_;
r_c = lambda*r(q);
while sqrt(sum((rd-r_c).^2)) > 0.001
    u = pinv(J(q))*(rd-r_c)+(eye(6)-pinv(J(q))*J(q))*(q_default-q);
    q = q+dt*u;
    r_c = lambda*r(q);
end
q0 = q;
% h = figure(4);hold on;
% c = [0.5,0.5,0.5];
% L = [1 1 1];
% plot_arm_3link (h,q0(1:3), L, c) ;axis tight;
% plot_arm_3link (h,q0(4:6), L, c) ;axis tight;
end

function [centres,s2,X] = read_to_matlab(filename,dim_x,dim_k,dim_r,dim_b)
fid = fopen(filename,'r');
formatSpec = '%f, ';
data = fscanf(fid,formatSpec);
c = 1;
for i = 1:dim_x
    for j = 1:dim_b
        centres(i,j) = data(c);
        c = c+1;
    end
end
s2 = data(c);
c = c+1;
for k = 1:dim_k
    for i = 1:dim_r-dim_k
        for j = 1:dim_b
            X{k}(i,j) = data(c);
            c = c+1;
        end
    end
end
fclose(fid);
end

function Lambda = predict_lambda (q, model, Iu)
    Rn      = Iu ;                                  % Initial rotation matrix
    Lambda  = zeros(model.dim_k, model.dim_r) ;     % Initial selection matrix
        
    for k = 1:model.dim_k              
        theta       = [pi/2 * ones(1,k-1) ,  (model.w{k} * model.phi(q) )' ] ;   
        alpha       = get_unit_vector_from_matrix(theta) ;                      % the kth alpha_0       
        Lambda(k,:) = alpha * Rn ;                                  % rotate alpha_0 to get the kth constraint
        Rn          = get_rotation_matrix (theta, Rn, model, k) ;   % update rotation matrix for (k+1)
    end                       
end
function R = get_rotation_matrix(theta_k, current_Rn, search, alpha_id)
    R = eye(search.dim_r) ;        
    for d = alpha_id : search.dim_t            
        R = R * make_givens_matrix(search.dim_r, d, d+1, theta_k(d) )' ;
    end        
    R = R * current_Rn ;
end
function G = make_givens_matrix(dim, i, j, theta)
    G      = eye(dim) ;
    G(i,i) = cos(theta) ;
    G(j,j) = cos(theta) ;
    G(i,j) =-sin(theta) ;
    G(j,i) = sin(theta) ;
end
function alpha = get_unit_vector_from_matrix (Theta)           
    [dim_n dim_t]   = size(Theta) ;      
    alpha           = zeros(dim_n,dim_t+1) ;     
    alpha(:,1)      = cos(Theta(:,1)) ;            
    for i =2:dim_t 
        alpha(:,i) = cos(Theta(:,i)) ;  
        
        for k = 1:i-1
            alpha(:,i) = alpha(:,i) .* sin(Theta(:,k)) ;        
        end                   
    end    
    alpha(:,dim_t+1)    = ones(dim_n,1) ;    
    for k = 1:dim_t            
        alpha(:,dim_t+1) = alpha(:,dim_t+1) .* sin(Theta(:,k)) ;    
    end        
end
