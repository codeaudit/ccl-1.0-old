function demo_resolve_nullspace_resolution
rng('default');
settings.dim_x = [3];
settings.dim_u = [3];
settings.dim_r = [2];
settings.dt = [0.02];
settings.N = [20];
[X,Y,Un] = data_generation(settings);                                          % define the number of radial basis functions
%% Learn_ncl
model.num_basis = 20;
model.c     = generate_kmeans_centres (X, model.num_basis) ;      % generate a grid of basis functions
model.s2    = real(mean(mean(sqrt(distances(model.c, model.c))))^2) ;   % set the variance as the mean distance between centres
model.phi   = @(x)phi_gaussian_rbf ( x, model.c, model.s2 ); % normalised Gaussian rbfs
tic;
model       = learn_ncl (X, Y, model) ;   % learn the model
t = toc;
s = sprintf('Finishing learning in %f ms',1000*t);
disp(s);
f_ncl       = @(x) predict_ncl ( model, x ) ; % set up an inference function
NSp = f_ncl (X) ;
%% evaluate
[NUPEtr,vNUPEtr_,NUPEtr_] = get_nupe(Un, NSp);
%% Native approach
BX = model.phi(X);
model_ = model;
model_.w = [];
model_ = learn_model_dir ( model_, BX, Y);
f_naive = @(x)predict_linear(x,model_);
%% Ground truth policy
lambda = [1 0 0;0 1 0];
robot.L = [1,1,1];
J  = @(q)J_planar_3_link_arm(q,robot);
f_T = @(x)((eye(3)-pinv(lambda*J(x))*(lambda*J(x)))*(1*[0.1745,-0.1745,0.1745]'-x));
%% reproduction
settings.N = 2;
filename = 'ncl_model.txt';
c_model = model;
c_model.w = [];
[c_model.c,c_model.s2,c_model.w] = read_to_matlab(filename,settings.dim_x,settings.dim_u,model.num_basis);
c_model.phi   = @(x)phi_gaussian_rbf ( x, c_model.c, c_model.s2 ); % normalised Gaussian rbfs
f_ncl_c       = @(x) predict_ncl ( c_model, x ) ; % set up an inference function
NSp_test = f_ncl_c (X) ;
for i = 1:settings.N
    q0= [0.1745*rand,0.1745*rand+1.57,0.1745*rand]';
    rd = [2*rand-1,2*rand]';
    title = 'Reproductions uisng True nullspace component';
    data_reproduction(settings,f_T,q0,rd,title);hold on;pause;
    title = 'Reproductions uisng Naive nullspace component';
    data_reproduction(settings,f_naive,q0,rd,title);pause;
    title = 'Reproductions uisng Learned nullspace component';
    data_reproduction(settings,f_ncl,q0,rd,title);pause;
    title = 'Reproductions uisng Learned nullspace component in C';
    data_reproduction(settings,f_ncl_c,q0,rd,title); %% from c code
    pause;
    clf;
end
close all;
end

function  [X,Y,Un] = data_generation(settings)
N = settings.N;
dt = settings.dt;
model.L = [1,1,1];
fk = @(q)r_planar_3_link_arm(q,model);
fk_all = @(q)r_planar_3_link_arm_all(q,model);
J  = @(q)J_planar_3_link_arm(q,model);
pi_ = @(q)(1*([0.1745,-0.1745,0.1745]'-q));
lambda = [1 0 0;0 1 0];
X = [];Y = [];Un = [];
for i = 1:N
    rd = [2*rand-1,2*rand]';
    q0  = [0.1745*rand,0.1745*rand+1.57,0.1745*rand]';
    k = 1;
    X_i = []; Y_i = []; Un_i = [];R_i = [];
    n = 1;
    q = q0;
    while n < 41
        r = lambda*fk(q);
        A = lambda*J(q);
        pinvA = pinv(A);
        N = eye(3) - pinvA*A;
        u_ts = pinvA*(k*(rd-r));
        u_ns = N*pi_(q);
        u =  u_ts+ u_ns;
        X_i = [X_i q];
        Y_i = [Y_i u];
        Un_i = [Un_i,u_ns];
        R_i(:,:,n) = fk_all(q);
        q = q + dt*u;
        n = n + 1;
    end
%     if i <= 3
%         fig = figure(1);
%         visualise_move_3link(fig,R_i,X_i,model.L,'Resolving nullspace resolutions',[0,500,600,600]);
%     end
    X = [X,X_i];
    Y = [Y,Y_i];
    Un = [Un,Un_i];
end
end

function  data_reproduction(settings,f_ncl,q0,rd,tt)
N = settings.N;
dt = settings.dt;
model.L = [1,1,1];
fk = @(q)r_planar_3_link_arm(q,model);
fk_all = @(q)r_planar_3_link_arm_all(q,model);
J  = @(q)J_planar_3_link_arm(q,model);
lambda = [1 0 0;0 1 0];
k = 1;
X_i = [];R_i = [];
n = 1;
q = q0;
r   = lambda*fk(q0);
while sqrt(sum((rd-r).^2)) > 0.1 && n < 50
    r   = lambda*fk(q0);
    A = lambda*J(q);
    pinvA = pinv(A);
    N = eye(3) - pinvA*A;
    u_ts = pinvA*(k*(rd-r));
    u_ns = f_ncl(q);
    u =  u_ts+ u_ns;
    X_i = [X_i q];
    R_i(:,:,n) = fk_all(q);
    q = q + dt*u;
    n = n + 1;
end
fig = figure(2);
visualise_move_3link(fig,R_i,X_i,model.L,tt,[600,500,600,600]);
end

function model = learn_model_dir ( model, BX, U )

HS  = eye(model.num_basis);
g   = BX * U' ;
H   = BX * BX';

% do eigen-decomposition for inversion
[V,D]   = eig( H + 1e-8*HS);
ev      = diag(D);
ind     = find( ev > 1e-8);
V1      = V(:,ind);
pinvH1  = V1 * diag(ev(ind).^-1)*V1';
model.w = (pinvH1 * g) ;
end

function [centres,s2,X] = read_to_matlab(filename,dim_x,dim_y,dim_b)
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
for i = 1:dim_y
    for j = 1:dim_b
        X(i,j) = data(c);
        c = c+1;
    end
end
fclose(fid);
end