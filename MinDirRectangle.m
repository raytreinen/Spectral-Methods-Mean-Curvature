% Jonas Haug, Rachel Jewell, Ray Treinen, May 2024
% 
% Compute capillary surfaces on the rectangle with Dirichlet data. 
% The rectangle is [aaa,bbb]x[ccc,ddd]
% The Dirichlet data is given by g
% kappa is set according to the physical problem.
%
% This function needs Chebfun installed to run: chebfun.org
% The dependencies on chebfun are restricted to the generation of the
% spectral differentiation matrices and plotting.

%% Physical parameters

aaa = -1;
bbb = 1;
ccc = -1;
ddd = 1;

g = @(x,y) 0.1*(sin(4*pi*x)).^2 + 0.1*(sin(4*pi*y));

%% Computational parameters
N = 55;

new_tol = 1e-14;
bvp_tol = 1e-10;
ep = 1e-8;
MM = 100;

%% Computational building blocks

x = chebpts(N,[aaa;bbb]);
y = chebpts(N,[ccc;ddd]);
[xx,yy] = meshgrid(x,y);
xx = xx(:);
yy = yy(:);

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

% Initial Guess
u0 = ones(size(xx));

b = find(xx == aaa | xx == bbb | yy == ccc | yy == ddd);
inside = find(xx ~= aaa & xx ~= bbb & yy ~= ccc & yy ~= ddd);


M = @(v) (Dxx*v) + (Dyy*v) + (Dxx*v).*(Dy*v).^2 - ...
    2*(Dx*v).*(Dy*v).*(Dxy*v) + (Dyy*v).*(Dx*v).^2;
Nu = zeros(size(u0));
Mu = M(u0);
Mui = Mu(inside);
Nu(inside) = Mui;
Nu(b) = u0(b)-g(xx(b),yy(b));

%% Solving the problem
bvp_res = 1;
count1 = 0;
while ((count1 < MM) && (bvp_res > bvp_tol))
    new_res2 = 1;
    count2 = 0;
    while ((count2 < MM) && (new_res2 > new_tol))

        M = @(v) (Dxx*v) + (Dyy*v) + (Dxx*v).*(Dy*v).^2 - ...
            2*(Dx*v).*(Dy*v).*(Dxy*v) + (Dyy*v).*(Dx*v).^2;

        F = @(v) Dxx + Dyy + ((Dy*v).^2).*Dxx + 2*(Dy*v).*(Dxx*v).*Dy + ...
            ((Dx*v).^2).*Dyy + 2*(Dx*v).*(Dyy*v).*Dx - 2*(Dx*v).*(Dxy*v).*Dx ...
            - 2*(Dx*v).*(Dxy*v).*Dy- 2*(Dx*v).*(Dy*v).*Dxy;

        Mu = M(u0);
        Mui = Mu(inside);
        Nu(inside) = Mui;
        Nu(b) = u0(b) - g(xx(b),yy(b));
        L = F(u0);
        for ii=1:length(b)
            z = zeros(1,length(u0));
            L(b(ii),:) = z;
            L(b(ii),b(ii)) = 1;
        end
        %size(L)
        %size(Nu)
        du = -L\Nu;
        new_res2 = norm(du)/(norm(u0) + ep);
        u0 = u0 + du;
        count2 = count2 + 1;

    end

    M = @(v) (Dxx*v) + (Dyy*v) + (Dxx*v).*(Dy*v).^2 - ...
        2*(Dx*v).*(Dy*v).*(Dxy*v) + (Dyy*v).*(Dx*v).^2;
    Mu = M(u0);
    Mui = Mu(inside);
    Nu(inside) = Mui;
    Nu(b) = u0(b) - g(xx(b),yy(b));

    bvp_res = norm(Nu)/(norm(u0) + ep);

    if (bvp_res > bvp_tol)

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
            2*(Dx*v).*(Dy*v).*(Dxy*v) + (Dyy*v).*(Dx*v).^2;
        Mu = M(u0);
        Mui = Mu(inside);
        Nu(inside) = Mui;
        Nu(b) = u0(b) - g(xx(b),yy(b));
    end
    count1 = count1 + 1;

end

%% Plotting

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
hold on
plot(Miloc(1),Miloc(2),'ok','MarkerSize',10)
plot(Maloc(1),Maloc(2),'.k','MarkerSize',30)
xlabel('X', 'FontWeight', 'bold')
ylabel('Y', 'FontWeight', 'bold')
fontsize("increase")

% ugrad = sqrt((Dx*u0).^2 + (Dy*u0).^2);
% ugrad = reshape(ugrad,N,N);
% ugrad = chebfun2(ugrad,[aaa bbb ccc ddd]);
% figure(5)
% plot(ugrad)
% title('Grad u')
% xlabel('X', 'FontWeight', 'bold')
% ylabel('Y', 'FontWeight', 'bold')
% zlabel('Grad U', 'FontWeight', 'bold')
% fontsize("increase")
% 
% [Mi, Miloc] = min2(ugrad);
% [Ma, Maloc] = max2(ugrad);
% 
% uy =  1 + (Dy*u0).^2;
% uy = reshape(uy,N,N);
% uy = chebfun2(uy,[aaa bbb ccc ddd]);
% figure(10)
% plot(uy)
% title('1 + u_y')
% xlabel('X', 'FontWeight', 'bold')
% ylabel('Y', 'FontWeight', 'bold')
% zlabel('U', 'FontWeight', 'bold')
% fontsize("increase")
% 
% [Mi, Miloc] = min2(uy);
% [Ma, Maloc] = max2(uy);
% 
% ux =  1 + (Dx*u0).^2;
% ux = reshape(ux,N,N);
% ux = chebfun2(ux,[aaa bbb ccc ddd]);
% figure(11)
% plot(ux)
% title('1 + u_x')
% xlabel('X', 'FontWeight', 'bold')
% ylabel('Y', 'FontWeight', 'bold')
% zlabel('U', 'FontWeight', 'bold')
% fontsize("increase")
% 
% [Mi, Miloc] = min2(ux);
% [Ma, Maloc] = max2(ux);
% 
% umix = -2*(Dy*u0).*(Dx*u0);
% umix = reshape(umix,N,N);
% umix = chebfun2(umix,[aaa bbb ccc ddd]);
% figure(12)
% plot(umix)
% title('-2*u_x*u_y')
% xlabel('X', 'FontWeight', 'bold')
% ylabel('Y', 'FontWeight', 'bold')
% zlabel('U', 'FontWeight', 'bold')
% fontsize("increase")
% 
% [Mi, Miloc] = min2(umix);
% [Ma, Maloc] = max2(umix);
% 
% uyyfactor =  2*((Dyy*u0).*(Dx*u0) - (Dxy*u0).*(Dy*u0));
% uyyfactor = reshape(uyyfactor,N,N);
% uyyfactor = chebfun2(uyyfactor,[aaa bbb ccc ddd]);
% figure(13)
% plot(uyyfactor)
% title('2*((Dyy*u0).*ux - (Dxy*u0).*uy)')
% xlabel('X', 'FontWeight', 'bold')
% ylabel('Y', 'FontWeight', 'bold')
% zlabel('U', 'FontWeight', 'bold')
% fontsize("increase")
% 
% [Mi, Miloc] = min2(uyyfactor);
% [Ma, Maloc] = max2(uyyfactor);
% 
% uxxfactor =  2*((Dxx*u0).*(Dy*u0) - (Dxy*u0).*(Dx*u0));
% uxxfactor = reshape(uxxfactor,N,N);
% uxxfactor = chebfun2(uxxfactor,[aaa bbb ccc ddd]);
% figure(14)
% plot(uxxfactor)
% title('2*((Dyy*u0).*ux - (Dxy*u0).*uy)')
% xlabel('X', 'FontWeight', 'bold')
% ylabel('Y', 'FontWeight', 'bold')
% zlabel('U', 'FontWeight', 'bold')
% fontsize("increase")
% 
% [Mi, Miloc] = min2(uxxfactor);
% [Ma, Maloc] = max2(uxxfactor);
