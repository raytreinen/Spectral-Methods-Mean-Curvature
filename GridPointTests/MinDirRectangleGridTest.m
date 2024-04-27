% Nvec = (20:2:70)';
Nvec = 55;
relerrp = zeros(length(Nvec),1);
relerr = zeros(length(Nvec),1);

% x between aaa and bbb
aaa = -1;
bbb = 1;
% y between ccc and ddd
ccc = -1;
ddd = 1;

new_tol = 1e-14;
ep = 1e-8;
MM = 100;

for i = 1:length(Nvec)
    N = Nvec(i)

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

    u0 = ones(size(xx));
    b = find(xx == aaa | xx == bbb | yy == ccc | yy == ddd);
    inside = find(xx ~= aaa & xx ~= bbb & yy ~= ccc & yy ~= ddd);
    g = @(x,y) 0.1*(sin(4*pi*x)).^2 + 0.1*(sin(4*pi*y));

    M = @(v) Dxx*v + Dyy*v + (Dxx*v).*(Dy*v).^2 - ...
        2*(Dx*v).*(Dy*v).*(Dxy*v) + (Dyy*v).*(Dx*v).^2;
    Nu = zeros(size(u0));
    Mu = M(u0);
    Mui = Mu(inside);
    Nu(inside) = Mui;
    Nu(b) = u0(b) - g(xx(b),yy(b));

    new_res2 = 1;
    count2 = 0;
    while ((count2 < MM) && (new_res2 > new_tol))

        M = @(v) Dxx*v + Dyy*v + (Dxx*v).*(Dy*v).^2 - ...
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

        du = -L\Nu;
        new_res2 = norm(du)/(norm(u0) + ep);
        u0 = u0 + du;
        count2 = count2 + 1;

    end

    M = @(v) Dxx*v + Dyy*v + (Dxx*v).*(Dy*v).^2 - ...
        2*(Dx*v).*(Dy*v).*(Dxy*v) + (Dyy*v).*(Dx*v).^2;
    Mu = M(u0);
    Mui = Mu(inside);
    Nu(inside) = Mui;
    Nu(b) = u0(b) - g(xx(b),yy(b));

    relerr(i) = norm(Nu)/(norm(u0) + ep);

    uu0 = reshape(u0, N, N);
    [xxc,yyc] = meshgrid(x,y);
    uu0 = chebfun2(uu0, [aaa bbb ccc ddd]);

    Np = N - 2;
    xp = chebpts(Np,[aaa;bbb]);
    yp = chebpts(Np,[ccc;ddd]);
    [xxp,yyp] = meshgrid(xp,yp);

    u0p = uu0(xxp,yyp);
    xxp = xxp(:);
    yyp = yyp(:);
    u0p = u0p(:);
    
    Ip = eye(Np);
    Dx1p = diffmat(Np,1,[aaa;bbb]);
    Dx2p = diffmat(Np,2,[aaa;bbb]);
    Dy1p = diffmat(Np,1,[ccc;ddd]);
    Dy2p = diffmat(Np,2,[ccc;ddd]);
    Dxp = kron(Dx1p,Ip);
    Dxxp = kron(Dx2p,Ip);
    Dyp = kron(Ip, Dy1p);
    Dyyp = kron(Ip, Dy2p);
    Dxyp = Dxp*Dyp;
    Dyxp = Dyp*Dxp;

    bp = find(xxp == aaa | xxp == bbb | yyp == ccc | yyp == ddd);
    insidep = find(xxp ~= aaa & xxp ~= bbb & yyp ~= ccc & yyp ~= ddd);

    Mp = @(v) Dxxp*v + Dyyp*v + (Dxxp*v).*(Dyp*v).^2 - ...
        2*(Dxp*v).*(Dyp*v).*(Dxyp*v) + (Dyyp*v).*(Dxp*v).^2;
    Mup = Mp(u0p);
    Muip = Mup(insidep);
    Nup = zeros(Np,1);
    Nup(insidep) = Muip;
    Nup(bp) = u0p(bp) - g(xxp(bp),yyp(bp));
    relerrp(i) = norm(Nup)/(norm(u0p) + ep);

    % Nui = norm(Nu(inside))
    % Nub = norm(Nu(b))
    % Nupi = norm(Nup(insidep))
    % Nupb = norm(Nup(bp))

    % figure(3)
    % plot(uu0)
    % xlabel('X', 'FontWeight', 'bold')
    % ylabel('Y', 'FontWeight', 'bold')
    % zlabel('U', 'FontWeight', 'bold')
    % fontsize("increase")
    % 
    % figure(5)
    % u0p = reshape(u0p, Np, Np);
    % u0p = chebfun2(u0p, [aaa bbb ccc ddd]);
    % plot(u0p)
    % xlabel('X', 'FontWeight', 'bold')
    % ylabel('Y', 'FontWeight', 'bold')
    % zlabel('U', 'FontWeight', 'bold')
    % fontsize("increase")

    % figure(4)
    % contour(xxc,yyc,uu0)
    % xlabel('X', 'FontWeight', 'bold')
    % ylabel('Y', 'FontWeight', 'bold')
    % fontsize("increase")
    % title(['Number of grid points is N = ',num2str(N)])
    % 
    % pause

    Nufig = reshape(Nu,N,N);
    Nufig = chebfun2(Nufig, [aaa bbb ccc ddd]);
    figure(12)
    plot(Nufig)
    xlabel('X', 'FontWeight', 'bold')
    ylabel('Y', 'FontWeight', 'bold')
    zlabel('U', 'FontWeight', 'bold')
    fontsize("increase")

    % Nupfig = reshape(Nup,Np,Np);
    % Nupfig = chebfun2(Nupfig, [aaa bbb ccc ddd]);
    % figure(13)
    % plot(Nupfig)
    % xlabel('X', 'FontWeight', 'bold')
    % ylabel('Y', 'FontWeight', 'bold')
    % zlabel('U', 'FontWeight', 'bold')
    % fontsize("increase")

    % pause
end

figure(7)
semilogy(Nvec,relerr,'g')% ,Nvec,relerrp,'k')
grid on
xlabel('Number of grid points', 'FontWeight', 'bold')
ylabel('Relative Error of N(u)', 'FontWeight', 'bold')
fontsize("increase")


%pause;

% uu0 = reshape(u0, N, N);
% [xx,yy] = meshgrid(x,y);
%
% uu0 = chebfun2(uu0, [aaa bbb ccc ddd]);
% figure(3)
% plot(uu0)
% xlabel('X', 'FontWeight', 'bold')
% ylabel('Y', 'FontWeight', 'bold')
% zlabel('U', 'FontWeight', 'bold')
% fontsize("increase")
% axis equal


