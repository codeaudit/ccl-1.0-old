function demo_resolve_nullspace_policy
rng('default');
settings.dim_x = [3];
settings.dim_u = [3];
settings.dim_r = [1];
settings.dt = [0.02];
settings.N = [20];
[X,Y,Pi] = data_generation(settings);                                          % define the number of radial basis functions
%% Learn_ncl
model.num_basis = 20;
model.c     = generate_kmeans_centres (X, model.num_basis) ;      % generate a grid of basis functions
model.s2    = real(mean(mean(sqrt(distances(model.c, model.c))))^2) ;   % set the variance as the mean distance between centres
model.W   = @(x)phi_gaussian_rbf( x, model.c, model.s2 );
model.phi = @(x)phi_linear( x );
model = learn_lwccl(X,Y,model); 
fp = @(x)predict_local_linear(x,model);
Fptr = fp(X);
%% evaluate
NUPEtr = get_nupe(Pi,Fptr);
%% Ground true nullspace policy
ft = @(q)(1*([0.1745,-0.1745,0.1745]'-q));
%% reproduction
settings.N = 2;
filename = 'learn_lwpi_model.txt';
c_model = model;
c_model.w = [];
[c_model.c,c_model.s2,c_model.w] = read_to_matlab(filename,settings.dim_x ,model.num_basis,settings.dim_u+1,settings.dim_u);
fp_c = @(x)predict_local_linear(x,c_model);
for i = 1:settings.N
    rad = pi*rand;
    lambda = [sin(rad),cos(rad),0];
    q0= [0.1745*rand,0.1745*rand+1.57,0.1745*rand]';
    title = 'Reproductions uisng True policy';
    data_reproduction(settings,ft,q0,lambda,title);hold on; pause;
    title = 'Reproductions uisng Learned policy';
    data_reproduction(settings,fp,q0,lambda,title); pause;
    title = 'Reproductions uisng Learned policy in C code';
    data_reproduction(settings,fp_c,q0,lambda,title); %% from c code
     pause;
    clf;
end
close all;
end

function  [X,Y,Pi] = data_generation(settings)
N = settings.N;
dt = settings.dt;
model.L = [1,1,1];
fk = @(q)r_planar_3_link_arm(q,model);
fk_all = @(q)r_planar_3_link_arm_all(q,model);
J  = @(q)J_planar_3_link_arm(q,model);
pi_ = @(q)(1*([0.1745,-0.1745,0.1745]'-q));
X = [];Y = [];Pi = [];
for i = 1:N
    rad = pi*rand;
    lambda = [sin(rad),cos(rad),0];
    q0  = [0.1745*rand,0.1745*rand+1.57,0.1745*rand]';
    k = 1;
    X_i = []; Y_i = []; Pi_i = [];R_i = [];
    n = 1;
    q = q0;
    while n < 41
        A = lambda*J(q);
        pinvA = pinv(A);
        N = eye(3) - pinvA*A;
        u_ts = zeros(3,1);
        u_ns = N*pi_(q);
        u =  u_ts+ u_ns;
        X_i = [X_i q];
        Y_i = [Y_i u];
        Pi_i = [Pi_i,pi_(q)];
        R_i(:,:,n) = fk_all(q);
        q = q + dt*u;
        n = n + 1;
    end
%     if i <= 3
%         fig = figure(1);
%         visualise_move_3link(fig,R_i,X_i,model.L,'Resolving nullspace policy -- a Wiping task',[0,500,600,600]);
%     end
    X = [X,X_i(:,2:end)];
    Y = [Y,Y_i(:,2:end)];
    Pi = [Pi,Pi_i(:,2:end)];
end
end

function  data_reproduction(settings,pi_hat,q0,lambda,tt)
dt = settings.dt;
model.L = [1,1,1];
fk_all = @(q)r_planar_3_link_arm_all(q,model);
J  = @(q)J_planar_3_link_arm(q,model);
X_i = [];R_i = [];
n = 1;
q = q0;
R_i = fk_all(q0);
X_i = q;
while n < 41
    A = lambda*J(q);
    pinvA = pinv(A);
    N = eye(3) - pinvA*A;
    u_ts =zeros(3,1);
    u_ns = N*pi_hat(q);
    u =  u_ts+ u_ns;
    q = q + dt*u;
    X_i = [X_i q];
    n = n + 1;
    R_i(:,:,n) = fk_all(q);
end
fig = figure(2);
visualise_move_3link(fig,R_i,X_i,model.L,tt,[600,500,600,600]);
end
function [centres,s2,X] = read_to_matlab(filename,dim_x,dim_b,dim_phi,dim_y)
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
for k = 1:dim_b
    for i = 1:dim_phi
        for j = 1:dim_y
            X(i,j,k) = data(c);
            c = c+1;
            if j == dim_y
                fprintf(fid,'\n');
            end
        end
    end
end
fclose(fid);
end