N = 30;
aaa = -1;
bbb = 1;
% y between ccc and ddd
ccc = -2;
ddd = 2;
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
% H is actually 2H
H = 1;

gamma = pi/2;% + 0.1; % 0.035;
cosg = cos(gamma);

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
bu = find(aaa < xx & xx < bbb & yy == ddd);
bd = find(aaa < xx & xx < bbb & yy== ccc);
bl = find(xx == aaa & ccc < yy & yy < ddd);
br = find(xx == bbb & ccc < yy & yy < ddd);
bur = find(xx == bbb & yy == ddd);
bul = find(xx == aaa & yy == ddd);
bdr = find(xx == bbb & yy == ccc);
bdl = find(xx == aaa & yy == ccc);

inside = find(xx ~= aaa & xx ~= bbb & yy ~= ccc & yy ~= ddd);

Dxu = Dx*u0;
Dyu = Dy*u0;

Dxbu = Dxu(bu);
Dybu = Dyu(bu);

Dxbd = Dxu(bd);
Dybd = Dyu(bd);

Dxbl = Dxu(bl);
Dybl = Dyu(bl);

Dxbr = Dxu(br);
Dybr = Dyu(br);

Dxbur = Dxu(bur);
Dybur = Dyu(bur);

Dxbul = Dxu(bul);
Dybul = Dyu(bul);

Dxbdr = Dxu(bdr);
Dybdr = Dyu(bdr);

Dxbdl = Dxu(bdl);
Dybdl = Dyu(bdl);

K = 100;
lambda = linspace(0,1,K)';
for j = 1:K

    M = @(v) (Dxx*v) + (Dyy*v) + (Dxx*v).*(Dy*v).^2 - ...
        2*(Dx*v).*(Dy*v).*(Dxy*v) + (Dyy*v).*(Dx*v).^2 - ...
        ((1 + ((Dx*v).^2) + ((Dy*v).^2)).^(3/2).*lambda(j)*H);
    Nu = zeros(size(u0));
    Mu = M(u0);
    Mui = Mu(inside);
    Nu(inside) = Mui;
    Nu(bu) = Dybu - lambda(j)*cosg*sqrt(1 + Dxbu.^2 + Dybu.^2);
    Nu(bd) = -Dybd - lambda(j)*cosg*sqrt(1 + Dxbd.^2 + Dybd.^2);
    Nu(bl) = -Dxbl - lambda(j)*cosg*sqrt(1 + Dxbl.^2 + Dybl.^2);
    Nu(br) = Dxbr - lambda(j)*cosg*sqrt(1 + Dxbr.^2 + Dybr.^2);
    Nu(bur) = Dxbur + Dybur - lambda(j)*2*cosg*sqrt(1 + Dxbur.^2 + Dybur.^2);
    Nu(bul) = -Dxbul + Dybul - lambda(j)*2*cosg*sqrt(1 + Dxbul.^2 + Dybul.^2);
    Nu(bdr) = Dxbdr - Dybdr - lambda(j)*2*cosg*sqrt(1 + Dxbdr.^2 + Dybdr.^2);
    Nu(bdl) = -Dxbdl - Dybdl - lambda(j)*2*cosg*sqrt(1 + Dxbdl.^2 + Dybdl.^2);


    new_res = 1;
    count2 = 0;
    while ((count2 < MM) && (new_res > new_tol))
    
        M = @(v) (Dxx*v) + (Dyy*v) + (Dxx*v).*(Dy*v).^2 - ...
            2*(Dx*v).*(Dy*v).*(Dxy*v) + (Dyy*v).*(Dx*v).^2 - ...
            ((1 + ((Dx*v).^2) + ((Dy*v).^2)).^(3/2).*lambda(j)*H);
    
        F = @(v) Dxx + Dyy + ((Dy*v).^2).*Dxx + 2*(Dy*v).*(Dxx*v).*Dy + ...
            ((Dx*v).^2).*Dyy + 2*(Dx*v).*(Dyy*v).*Dx - 2*(Dx*v).*(Dxy*v).*Dx ...
            - 2*(Dx*v).*(Dxy*v).*Dy- 2*(Dx*v).*(Dy*v).*Dxy - ...
            3 * lambda(j)*H * (1 + ((Dx*v).^2) + ((Dy*v).^2)).^(1/2) .* ...
            (((Dx*v).*Dx) + ((Dy*v).*Dy));
    
        Mu = M(u0);
        Mui = Mu(inside);
        Nu(inside) = Mui;
        L = F(u0);
        
        Dxu = Dx*u0;
        Dyu = Dy*u0;

        Dxbu = Dxu(bu);
        Dybu = Dyu(bu);
        L(bu,:) = (-lambda(j)*cosg*Dxbu./sqrt(1 + Dxbu.^2 + Dybu.^2)).*Dx(bu,:) + (ones(length(bu),1) - lambda(j)*cosg*Dybu./sqrt(1 + Dxbu.^2 + Dybu.^2)).*Dy(bu,:);
        Nu(bu) = Dybu - lambda(j)*cosg*sqrt(1 + Dxbu.^2 + Dybu.^2);

        Dxbd = Dxu(bd);
        Dybd = Dyu(bd);
        L(bd,:) = (-lambda(j)*cosg*Dxbu./sqrt(1 + Dxbd.^2 + Dybd.^2)).*Dx(bd,:) + (-ones(length(bd),1) - lambda(j)*cosg*Dybd./sqrt(1 + Dxbd.^2 + Dybd.^2)).*Dy(bd,:);
        Nu(bd) = -Dybd - lambda(j)*cosg*sqrt(1 + Dxbd.^2 + Dybd.^2);

        Dxbl = Dxu(bl);
        Dybl = Dyu(bl);
        L(bl,:) = (-ones(length(bl),1)-lambda(j)*cosg*Dxbl./sqrt(1 + Dxbl.^2 + Dybl.^2)).*Dx(bl,:) + (-lambda(j)*cosg*Dybl./sqrt(1 + Dxbl.^2 + Dybl.^2)).*Dy(bl,:);
        Nu(bl) = -Dxbl - lambda(j)*cosg*sqrt(1 + Dxbl.^2 + Dybl.^2);

        Dxbr = Dxu(br);
        Dybr = Dyu(br);
        L(br,:) = (ones(length(br),1)-lambda(j)*cosg*Dxbr./sqrt(1 + Dxbr.^2 + Dybr.^2)).*Dx(br,:) + (-lambda(j)*cosg*Dybr./sqrt(1 + Dxbr.^2 + Dybr.^2)).*Dy(br,:);
        Nu(br) = Dxbr - lambda(j)*cosg*sqrt(1 + Dxbr.^2 + Dybr.^2);

        Dxbur = Dxu(bur);
        Dybur = Dyu(bur);
        L(bur,:) = (1- lambda(j)*2*cosg*Dxbur./sqrt(1 + Dxbur.^2 + Dybur.^2)).*Dx(bur,:) + (1 -lambda(j)*2*cosg*Dybur./sqrt(1 + Dxbur.^2 + Dybur.^2)).*Dy(bur,:);
        Nu(bur) = Dxbur + Dybur - lambda(j)*2*cosg*sqrt(1 + Dxbur.^2 + Dybur.^2);

        Dxbul = Dxu(bul);
        Dybul = Dyu(bul);
        L(bul,:) = (-1- lambda(j)*2*cosg*Dxbul./sqrt(1 + Dxbul.^2 + Dybul.^2)).*Dx(bul,:) + (1 -lambda(j)*2*cosg*Dybul./sqrt(1 + Dxbul.^2 + Dybul.^2)).*Dy(bul,:);
        Nu(bul) = -Dxbul + Dybul - lambda(j)*2*cosg*sqrt(1 + Dxbul.^2 + Dybul.^2);

        Dxbdr = Dxu(bdr);
        Dybdr = Dyu(bdr);
        L(bdr,:) = (1- lambda(j)*2*cosg*Dxbdr./sqrt(1 + Dxbdr.^2 + Dybdr.^2)).*Dx(bdr,:) + (-1 -lambda(j)*2*cosg*Dybdr./sqrt(1 + Dxbdr.^2 + Dybdr.^2)).*Dy(bdr,:);
        Nu(bdr) = Dxbdr - Dybdr - lambda(j)*2*cosg*sqrt(1 + Dxbdr.^2 + Dybdr.^2);

        Dxbdl = Dxu(bdl);
        Dybdl = Dyu(bdl);
        L(bdl,:) = (-1- lambda(j)*2*cosg*Dxbdl./sqrt(1 + Dxbdl.^2 + Dybdl.^2)).*Dx(bdl,:) + (-1 -lambda(j)*2*cosg*Dybdl./sqrt(1 + Dxbdl.^2 + Dybdl.^2)).*Dy(bdl,:);
        Nu(bdl) = -Dxbdl - Dybdl - lambda(j)*2*cosg*sqrt(1 + Dxbdl.^2 + Dybdl.^2);

        du = -L\Nu;
        new_res = norm(du) / (norm(u0)+ep);
        u0 = u0 + du;
        count2 = count2 + 1;
    
    end
    uu0 = reshape(u0, N, N);

    uu0 = chebfun2(uu0, [aaa bbb ccc ddd]);
    figure(7)
    plot(uu0)
    drawnow
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
axis equal

figure(4)
contour(xx,yy,uu0)
xlabel('X', 'FontWeight', 'bold')
ylabel('Y', 'FontWeight', 'bold')
fontsize("increase")
