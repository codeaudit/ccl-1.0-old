
function Yp = predict_local_linear(X,model)
dim_x = size(X,1);
W = model.W(X);
nc = size(W,1);
Phi = model.phi(X);
Yp = zeros(size(X));
for i = 1:nc
    Yp  = Yp + (Phi'*model.w(:,:,i).*repmat(W(i,:)',1,dim_x))';
end
