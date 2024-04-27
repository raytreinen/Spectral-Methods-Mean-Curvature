N = 55;
aaa = -1;
bbb = 1;
% y between ccc and ddd
ccc = -1;
ddd = 1;
x = chebpts(N,[aaa;bbb]);
% y = x;
y = chebpts(N,[ccc;ddd]);
[xx,yy] = meshgrid(x,y);
xx = xx(:);
yy = yy(:);
new_tol = 1e-14;
bvp_tol = 1e-10;
ep = 1e-8;
MM = 100;
%will start small, then increase 
kappa = 1;

I = eye(N);
Dx1 = diffmat(N,1,[aaa;bbb]);
Dx2 = diffmat(N,2,[aaa;bbb]);
Dy1 = diffmat(N,1,[ccc;ddd]);
Dy2 = diffmat(N,2,[ccc;ddd]);
Dx = kron(Dx1,I);
Dxx = kron(Dx2,I);
Dy = kron(I, Dy1);
Dyy = kron(I, Dy2);
Dxy = Dx * Dy;
Dyx = Dy * Dx;

u0 = zeros(size(xx));
b = find(xx == aaa | xx == bbb | yy == ccc | yy == ddd);
inside = find(xx ~= aaa & xx ~= bbb & yy ~= ccc & yy ~= ddd);
g = @(x,y) 0.1*sin(2*pi*x).^2 + 0.2*sin(1*pi*y) + 1;

M = @(v) (Dxx*v) + (Dyy*v) + (Dxx*v).*(Dy*v).^2 - ...
    2*(Dx*v).*(Dy*v).*(Dxy*v) + (Dyy*v).*(Dx*v).^2 - ...
    ((1 + ((Dx*v).^2) + ((Dy*v).^2)).^(3/2).*kappa.*v);
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
            ((1 + ((Dx*v).^2) + ((Dy*v).^2)).^(3/2).*kappa.*v);
    
        F = @(v) Dxx + Dyy + ((Dy*v).^2).*Dxx + 2*(Dy*v).*(Dxx*v).*Dy + ...
            ((Dx*v).^2).*Dyy + 2*(Dx*v).*(Dyy*v).*Dx - 2*(Dx*v).*(Dxy*v).*Dx ...
            - 2*(Dx*v).*(Dxy*v).*Dy- 2*(Dx*v).*(Dy*v).*Dxy - ...
            3 * kappa.*v .* (1 + ((Dx*v).^2) + ((Dy*v).^2)).^(1/2) .* ...
            (((Dx*v).*Dx) + ((Dy*v).*Dy)) - ...
            (1 + ((Dx*v).^2) + ((Dy*v).^2)).^(3/2).*kappa .* eye(N^2);
    
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
            ((1 + ((Dx*v).^2) + ((Dy*v).^2)).^(3/2).*kappa.*v);
    Mu = M(u0);
    Mui = Mu(inside);
    Nu(inside) = Mui;
    Nu(b) = u0(b)-g(xx(b),yy(b));
    
    bvp_res = (norm(Nu)) / (norm(u0)+ep);
    
    if (bvp_res > new_tol)

        uu0 = reshape(u0, N, N);
        uu0 = chebfun2(uu0, [aaa bbb ccc ddd]);

        N = N + 4;
        I = eye(N);
        Dx1 = diffmat(N,1,[aaa;bbb]);
        Dx2 = diffmat(N,2,[aaa;bbb]);
        Dy1 = diffmat(N,1,[ccc;ddd]);
        Dy2 = diffmat(N,2,[ccc;ddd]);
        Dx = kron(Dx1,I);
        Dxx = kron(Dx2,I);
        Dy = kron(I, Dy1);
        Dyy = kron(I, Dy2);

        Dxy = Dx * Dy;
        Dyx = Dy * Dx;
        
        x = chebpts(N,[aaa;bbb]);
        y = chebpts(N,[ccc;ddd]);
        [xx,yy] = meshgrid(x,y);

        u0 = uu0(xx,yy);
        xx = xx(:);
        yy = yy(:);
        u0 = u0(:);

        b = find(xx == aaa | xx == bbb | yy == ccc | yy == ddd);
        inside = find(xx ~= aaa & xx ~= bbb & yy ~= ccc & yy ~= ddd);

        M = @(v) (Dxx*v) + (Dyy*v) + (Dxx*v).*(Dy*v).^2 - ...
            2*(Dx*v).*(Dy*v).*(Dxy*v) + (Dyy*v).*(Dx*v).^2 - ...
            ((1 + ((Dx*v).^2) + ((Dy*v).^2)).^(3/2).*kappa.*v);
        Mu = M(u0);
        Mui = Mu(inside);
        Nu(inside) = Mui;
        Nu(b) = u0(b) - g(xx(b),yy(b));
    end
    count1 = count1 + 1;

end

uu0 = reshape(u0, N, N);
[xx,yy] = meshgrid(x,y);

uu0 = chebfun2(uu0, [aaa bbb ccc ddd]);
figure(3)
plot(uu0)
xlabel('X', 'FontWeight', 'bold')
ylabel('Y', 'FontWeight', 'bold')
zlabel('U', 'FontWeight', 'bold')
fontsize("increase")
% axis equal

figure(4)
contour(xx,yy,uu0)
xlabel('X', 'FontWeight', 'bold')
ylabel('Y', 'FontWeight', 'bold')
fontsize("increase")

hold on
plot([aaa;bbb], [ccc;ccc], 'k')
plot([aaa;bbb], [ddd;ddd], 'k')
plot([aaa;aaa], [ccc;ddd], 'k')
plot([bbb;bbb], [ccc;ddd], 'k')