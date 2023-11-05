N = 45;
x = chebpts(N);
y = x;
[xx,yy] = meshgrid(x,y);
xx = xx(:);
yy = yy(:);
new_tol = 1e-14;
bvp_tol = 1e-10;
ep = 1e-8;
MM = 100;
%will start small, then increase 
lambda = 0.5;

I = eye(N);
D1 = diffmat(N,1);
D2 = diffmat(N,2);
Dx = kron(D1,I);
Dxx = kron(D2,I);
Dy = kron(I, D1);
Dyy = kron(I, D2);
Dxy = Dx * Dy;
Dyx = Dy * Dx;

u0 = zeros(size(xx));
b = find(abs(xx)==1 | abs(yy)==1);
inside = find(abs(xx)~=1 & abs(yy)~=1);
g = @(x,y) 0.1*sin(2*pi*x).^2;
g = @(x,y) 0.1*sin(2*pi*x).^2 + 0.1*sin(4*pi*y);

M = @(v) (Dxx*v) + (Dyy*v) + (Dxx*v).*(Dy*v).^2 - ...
    2*(Dx*v).*(Dy*v).*(Dxy*v) + (Dyy*v).*(Dx*v).^2 - ...
    ((1 + ((Dx*v).^2) + ((Dy*v).^2)).^(3/2).*lambda);
Nu = zeros(size(u0));
Mu = M(u0);
Mui = Mu(inside);
Nu(inside) = Mui;
Nu(b) = u0(b)-g(xx(b),yy(b));

bvp_res = 1;
count1 = 0;
while ((count1 < MM) && (bvp_res > bvp_tol))
    new_res = 1;
    count2 = 0;
    while ((count2 < MM) && (new_res > new_tol))
    
        M = @(v) (Dxx*v) + (Dyy*v) + (Dxx*v).*(Dy*v).^2 - ...
            2*(Dx*v).*(Dy*v).*(Dxy*v) + (Dyy*v).*(Dx*v).^2 - ...
            ((1 + ((Dx*v).^2) + ((Dy*v).^2)).^(3/2).*lambda);
    
        F = @(v) Dxx + Dyy + ((Dy*v).^2).*Dxx + 2*(Dy*v).*(Dxx*v).*Dy + ...
            ((Dx*v).^2).*Dyy + 2*(Dx*v).*(Dyy*v).*Dx - 2*(Dx*v).*(Dxy*v).*Dx ...
            - 2*(Dx*v).*(Dxy*v).*Dy- 2*(Dx*v).*(Dy*v).*Dxy - ...
            3 * lambda * (1 + ((Dx*v).^2) + ((Dy*v).^2)).^(1/2) .* ...
            (((Dx*v).*Dx) + ((Dy*v).*Dy));
    
        Mu = M(u0);
        Mui = Mu(inside);
        Nu(inside) = Mui;
        Nu(b) = u0(b)-g(xx(b),yy(b));
        L = F(u0);
        for ii=1:length(b)
            z = zeros(1,length(u0));
            L(b(ii),:) = z;
            L(b(ii),b(ii)) = 1;
        end
        
        du = -L\Nu;
        new_res = (norm(du)) / (norm(u0)+ep);
        u0 = u0+du;
        count2 = count2 + 1;
    
    end
    
    M = @(v) (Dxx*v) + (Dyy*v) + (Dxx*v).*(Dy*v).^2 - ...
            2*(Dx*v).*(Dy*v).*(Dxy*v) + (Dyy*v).*(Dx*v).^2 - ...
            ((1 + ((Dx*v).^2) + ((Dy*v).^2)).^(3/2).*lambda);
    Mu = M(u0);
    Mui = Mu(inside);
    Nu(inside) = Mui;
    Nu(b) = u0(b)-g(xx(b),yy(b));
    
    bvp_res = (norm(Nu)) / (norm(u0)+ep);
    
    if (bvp_res > new_tol)

        uu0 = reshape(u0, N, N);
        uu0 = chebfun2(uu0);

        N = N + 4;
        I = eye(N);
        D1 = diffmat(N,1);
        D2 = diffmat(N,2);
        Dx = kron(D1,I);
        Dxx = kron(D2,I);
        Dy = kron(I, D1);
        Dyy = kron(I, D2);
        Dxy = Dx * Dy;
        Dyx = Dy * Dx;

       
        
        x = chebpts(N);
        y = x;
        [xx,yy] = meshgrid(x,y);

        u0 = uu0(xx,yy);
        xx = xx(:);
        yy = yy(:);
        u0 = u0(:);

        b = find(abs(xx)==1 | abs(yy)==1);
        inside = find(abs(xx)~=1 & abs(yy)~=1);
        
        M = @(v) (Dxx*v) + (Dyy*v) + (Dxx*v).*(Dy*v).^2 - ...
            2*(Dx*v).*(Dy*v).*(Dxy*v) + (Dyy*v).*(Dx*v).^2 - ...
            ((1 + ((Dx*v).^2) + ((Dy*v).^2)).^(3/2).*lambda);
        Mu = M(u0);
        Mui = Mu(inside);
        Nu(inside) = Mui;
        Nu(b) = u0(b)-g(xx(b),yy(b));
    end
    count1 = count1 + 1;

end

uu0 = reshape(u0, N, N);
[xx,yy] = meshgrid(x,y);

uu0 = chebfun2(uu0);
plot(uu0)