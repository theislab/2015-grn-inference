function [psi,E] = diffusion_maps_nn(data, no_dims, nn, t, sigma)

% data    = high dimensional matrix data to be analysied 
% nn      = number of neighbours to be considered in the analysis, default should be number of cells in the data size(data,1)
% no_dims = number of mapping dimensions, defualt should be where there is a gap between the eigenvalues E
% t       = density normalisation factor, default should be 1
% sigma   = diffusion scale parameter of the Gauusian kernel
% psi     = mapped coordinates of each cell

KT = 2*sigma^2;
n  = size(data,1);
W  = zeros(n,n);

d2 = zeros(n,n);
R  = zeros(n,1);

for i = 1:n
    for j = 1:n
        d2(i,j)=(dot(data(i,:)-data(j,:),data(i,:)-data(j,:)));
    end
    [sortd2,id] = sort(d2(i,:));
    R(i)        = sortd2(nn); 
end

for i = 1:n
    for j = 1:n
        if d2(i,j)<R(i)
            W(i,j) = exp(-d2(i,j)/KT);
        end
    end
    W(i,i) = 0;
end

D = zeros(n,n);
for i = 1:n
    D(i,i) = sum(W(:,i));
end
 
q = zeros(n,n);
for i = 1:n
    for j = 1:n
        q(i,j) = (D(i,i)*D(j,j))^(t);
    end
end

H  = W./q;
D_ = zeros(n,n);

for i = 1:n
    D_(i,i) = sum(H(:,i));
end

Hp = D_^(-1)*H;

[psi_nsort,En] = eig(Hp);
[E, indE]      = sort(real(diag(En)),'descend'); 
psi            = psi_nsort(:,indE); 
psi            = psi(:,2:no_dims+1);
E              = E(2:20);
    