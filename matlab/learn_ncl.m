%{
learn the nullspace component of the input data X

input
    X: observed states
    Y: observed actions of the form Y = A(X)'B(X) + N(X) F(X) where N and F 
       are consistent across the dataset X 

output
    model: model that prodicts the nullspace component N(X)F(X) 
%}
function model = learn_ncl(X, Y, model)
    
    % calculate the basis functions of the training data    
    BX      = model.phi(X) ;
    
    % learn an initial model using a simple parametric model
    model   = learn_model_dir ( model, BX, Y ) ;  
%     model   = learn_lwccl ( X, Y, model ) ;  
%     f_ncl   = @(x) predict_ncl ( model, x ) ;
%     y_hat    = f_ncl (X) ;
    % learn the model by minimising the proposed step-1 method        
    obj     = @(W) obj_ncl_reg ( model, W, BX, Y);   % setup the learning objective function 
%     obj     = @(B) obj_ncl_reg_lw ( model, B, X, Y);   % setup the learning objective function  
%     options = optimset( 'Jacobian','on', 'Display', 'notify',...  % options for the optimisation    
%                         'MaxFunEvals',1e9, 'MaxIter', 1000,...
%                         'TolFun',1e-9, 'TolX',1e-9,...
%                         'Algorithm', 'levenberg-marquardt');
%     model.w = lsqnonlin(obj, model.w, [], [], options );  % use the non-linear optimiser to solve obj_ncl_reg    

%% test
    options.MaxIter = 1000 ;
    options.TolFun  = 1e-9 ;
    options.TolX    = 1e-9 ;   
    model.w   = solve_lm_ncl (obj, model.w, options );
    model.w   = reshape(model.w,size(Y,1),size(BX,1));


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
    model.w = (pinvH1 * g)' ;    
end

function model = learn_lwccl(X,Y,model)

[dimY N] = size(Y);

% find normalised Y
r = sum(Y.^2,1).^0.5;
YN = Y./repmat(r,dimY,1);

% find feature vectors
Phi      = model.phi(X);
dimPhi   = size(Phi,1); % get feature dimensionality

% find weights
W        = model.W(X);
Nc       = size(W,1);   % get no. centres

% train each local model
for nc=1:Nc
	WPhi=repmat(W(nc,:),dimPhi,1).*Phi;

	% construct Jacobian
	YPhit = Y*WPhi';
	g = YPhit(:);

	% construct Hessian
	H = zeros(dimY*dimPhi);
	for n=1:N
	YNPhit = YN(:,n)*Phi(:,n)';
	v(:,n) = YNPhit(:);
	H = H + W(nc,n)*v(:,n)*v(:,n)';
	end

	% do eigendecomposition for inversion
	%[V,D] = eig(H+1e-6*eye(size(H)));
	[V,D] = eig(H);
	ev = diag(D);
	ind = find(ev>1e-6);
	V1=V(:,ind);
	pinvH1 = V1*diag(ev(ind).^-1)*V1';
	model.w(:,:,nc)=reshape(pinvH1*g,dimY,dimPhi)';
end
end
