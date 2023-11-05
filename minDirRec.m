N = 30;
x = chebpts(N);
y = x;
[xx,yy] = meshgrid(x,y);
xx = xx(:);
yy = yy(:);
new_tol = 1e-14;
bvp_tol = 1e-10;
ep = 1e-8;
MM = 100;

I = eye(N);
D1 = diffmat(N,1);
D2 = diffmat(N,2);
Dx = kron(D1,I);
Dxx = kron(D2,I);
Dy = kron(I, D1);
Dyy = kron(I, D2);
Dxy = Dx * Dy;
Dyx = Dy * Dx;

u0 = ones(size(xx));
b = find(abs(xx)==1 | abs(yy)==1);
inside = find(abs(xx)~=1 & abs(yy)~=1);
g = @(x,y) 0.1*(sin(4*pi*x)).^2 + 0.1*(sin(4*pi*y));

M = @(v) (Dxx*v) + (Dyy*v) + (Dxx*v).*(Dy*v).^2 - ...
    2*(Dx*v).*(Dy*v).*(Dxy*v) + (Dyy*v).*(Dx*v).^2;
Nu = zeros(size(u0));
Mu = M(u0);
Mui = Mu(inside);
Nu(inside) = Mui;
Nu(b) = u0(b)-g(xx(b),yy(b));

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
        Nu(b) = u0(b)-g(xx(b),yy(b));
        L = F(u0);
        for ii=1:length(b)
            z = zeros(1,length(u0));
            L(b(ii),:) = z;
            L(b(ii),b(ii)) = 1;
        end
        %size(L)
        %size(Nu)
        du = -L\Nu;
        new_res2 = (norm(du)) / (norm(u0)+ep);
        u0 = u0+du;
        count2 = count2 + 1;

    end

    M = @(v) (Dxx*v) + (Dyy*v) + (Dxx*v).*(Dy*v).^2 - ...
        2*(Dx*v).*(Dy*v).*(Dxy*v) + (Dyy*v).*(Dx*v).^2;
    Mu = M(u0);
    Mui = Mu(inside);
    Nu(inside) = Mui;
    Nu(b) = u0(b)-g(xx(b),yy(b));

    %     new_res2
    %         uu0 = reshape(u0, N, N);
    %         uu0 = chebfun2(uu0);
    %
    %         Nb = N + 4;
    %         Ib = eye(Nb);
    %         D1b = diffmat(Nb,1);
    %         D2b = diffmat(Nb,2);
    %         Dxb = kron(D1b,Ib);
    %         Dxxb = kron(D2b,Ib);
    %         Dyb = kron(Ib, D1b);
    %         Dyyb = kron(Ib, D2b);
    %         Dxyb = Dxb * Dyb;
    %         Dyxb = Dyb * Dxb;
    %
    %         xb = chebpts(Nb);
    %         yb = xb;
    %         [xxb,yyb] = meshgrid(xb,yb);
    %
    %         u0b = uu0(xxb,yyb);
    %         xxb = xxb(:);
    %         yyb = yyb(:);
    %         u0b = u0b(:);
    %
    %         bb = find(abs(xxb)==1 | abs(yyb)==1);
    %         insideb = find(abs(xxb)~=1 & abs(yyb)~=1);
    % %         g = xx(b).^2 + yy(b).^3;
    % %         size(u0b)
    % %         size(Dxxb)
    %         Mb = @(v) (Dxxb*v) + (Dyyb*v) + (Dxxb*v).*(Dyb*v).^2 - ...
    %         2*(Dxb*v).*(Dyb*v).*(Dxyb*v) + (Dyyb*v).*(Dxb*v).^2;
    %         Mub = Mb(u0b);
    %         Muib = Mub(insideb);
    %         Nub(insideb) = Muib;
    %         Nub(bb) = u0b(bb)-g(xxb(bb),yyb(bb));
    %
    %
    %     bvp_res = (norm(Nub)) / (norm(u0b)+ep)
    bvp_res = (norm(Nu)) / (norm(u0)+ep)
    %pause;
    if (bvp_res > bvp_tol)

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
        %         g = xx(b).^2 + yy(b).^3;

        M = @(v) (Dxx*v) + (Dyy*v) + (Dxx*v).*(Dy*v).^2 - ...
            2*(Dx*v).*(Dy*v).*(Dxy*v) + (Dyy*v).*(Dx*v).^2;
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