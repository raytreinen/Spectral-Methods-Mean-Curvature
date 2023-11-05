N = 40;
x = chebpts(N);
y = x;
[xx,yy] = meshgrid(x,y);
xx = xx(:);
yy = yy(:);
new_tol = 1e-12;
bvp_tol = 1e-10;
ep = 1e-8;
MM = 100;
kappa = 0.001;
kappa = 1;

gamma = pi/4 + 0.035;
cosg = cos(gamma);


I = eye(N);
D1 = diffmat(N,1);
D2 = diffmat(N,2);
Dx = kron(D1,I);
Dxx = kron(D2,I);
Dy = kron(I, D1);
Dyy = kron(I, D2);
Dxy = Dx * Dy;
Dyx = Dy * Dx;

% u0 = ones(size(xx));
[xxx,yyy] = meshgrid(x,y);
u0 = 2*sqrt(2) + 1/(2*sqrt(2)) - sqrt(8 - xxx.^2 - yyy.^2);
% figure(3)
% surf(xxx,yyy,u0)
% pause
u0 = u0(:);

% b = find(abs(xx)==1 | abs(yy)==1);
% bu = find(abs(xx)<1 | yy==1);
% bd = find(abs(xx)<1 | yy==-1);
% bl = find((xx)==-1 | abs(yy)<1);
% br = find((xx)==1 | abs(yy)<1);
% bur = find(xx==1 | yy==1);
% bul = find(xx==-1 | yy==1);
% bdr = find(xx==1 | yy==-1);
% bdl = find(xx==-1 | yy==-1);

b = find(abs(xx)==1 | abs(yy)==1);
bu = find(abs(xx)<1 & yy==1);
bd = find(abs(xx)<1 & yy==-1);
bl = find((xx)==-1 & abs(yy)<1);
br = find((xx)==1 & abs(yy)<1);
bur = find(xx==1 & yy==1);
bul = find(xx==-1 & yy==1);
bdr = find(xx==1 & yy==-1);
bdl = find(xx==-1 & yy==-1);

inside = find(abs(xx)~=1 & abs(yy)~=1);

M = @(v) (Dxx*v) + (Dyy*v) + (Dxx*v).*(Dy*v).^2 - ...
    2*(Dx*v).*(Dy*v).*(Dxy*v) + (Dyy*v).*(Dx*v).^2 - ...
    ((1 + ((Dx*v).^2) + ((Dy*v).^2)).^(3/2).*kappa.*v);
Nu = zeros(size(u0));
Mu = M(u0);
Mui = Mu(inside);
Nu(inside) = Mui;

Dxu = Dx*u0;
Dyu = Dy*u0;

Dxbu = Dxu(bu);
Dybu = Dyu(bu);
Nu(bu) = Dybu - cosg*sqrt(1 + Dxbu.^2 + Dybu.^2);

Dxbd = Dxu(bd);
Dybd = Dyu(bd);
Nu(bd) = -Dybd - cosg*sqrt(1 + Dxbd.^2 + Dybd.^2);

Dxbl = Dxu(bl);
Dybl = Dyu(bl);
Nu(bl) = -Dxbl - cosg*sqrt(1 + Dxbl.^2 + Dybl.^2);

Dxbr = Dxu(br);
Dybr = Dyu(br);
Nu(br) = Dxbr - cosg*sqrt(1 + Dxbr.^2 + Dybr.^2);

Dxbur = Dxu(bur);
Dybur = Dyu(bur);
%Nu(bur) = Dxbur + Dybur - 2*sqrt(2)*cosg*sqrt(1 + Dxbur.^2 + Dybur.^2);
Nu(bur) = Dxbur + Dybur - 2*cosg*sqrt(1 + Dxbur.^2 + Dybur.^2);

Dxbul = Dxu(bul);
Dybul = Dyu(bul);
%Nu(bul) = -Dxbul + Dybul - 2*sqrt(2)*cosg*sqrt(1 + Dxbul.^2 + Dybul.^2);
Nu(bul) = -Dxbul + Dybul - 2*cosg*sqrt(1 + Dxbul.^2 + Dybul.^2);

Dxbdr = Dxu(bdr);
Dybdr = Dyu(bdr);
%Nu(bdr) = Dxbdr - Dybdr - 2*sqrt(2)*cosg*sqrt(1 + Dxbdr.^2 + Dybdr.^2);
Nu(bdr) = Dxbdr - Dybdr - 2*cosg*sqrt(1 + Dxbdr.^2 + Dybdr.^2);

Dxbdl = Dxu(bdl);
Dybdl = Dyu(bdl);
%Nu(bdl) = Dxbdl + Dybdl + 2*sqrt(2)*cosg*sqrt(1 + Dxbdl.^2 + Dybdl.^2);
Nu(bdl) = -Dxbdl - Dybdl - 2*cosg*sqrt(1 + Dxbdl.^2 + Dybdl.^2);

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
        L = F(u0);
                 
        Dxu = Dx*u0;
        Dyu = Dy*u0;

        Dxbu = Dxu(bu);
        Dybu = Dyu(bu);
        % L(bu,:) = cosg*Dxbu.*Dx(bu,:) + (cosg*Dybu - sqrt(1 + Dxbu.^2 + Dybu.^2)).*Dy(bu,:);
        L(bu,:) = (-cosg*Dxbu./sqrt(1 + Dxbu.^2 + Dybu.^2)).*Dx(bu,:) + (ones(length(bu),1)-cosg*Dybu./sqrt(1 + Dxbu.^2 + Dybu.^2)).*Dy(bu,:);
        Nu(bu) = Dybu - cosg*sqrt(1 + Dxbu.^2 + Dybu.^2);

        Dxbd = Dxu(bd);
        Dybd = Dyu(bd);
        % L(bd,:) = cosg*Dxbd.*Dx(bd,:) + (cosg*Dybd + sqrt(1 + Dxbd.^2 + Dybd.^2)).*Dy(bd,:);
        L(bd,:) = (-cosg*Dxbu./sqrt(1 + Dxbd.^2 + Dybd.^2)).*Dx(bd,:) + (-ones(length(bd),1)-cosg*Dybd./sqrt(1 + Dxbd.^2 + Dybd.^2)).*Dy(bd,:);
        Nu(bd) = -Dybd - cosg*sqrt(1 + Dxbd.^2 + Dybd.^2);

        Dxbl = Dxu(bl);
        Dybl = Dyu(bl);
        % L(bl,:) = (cosg*Dxbl + sqrt(1 + Dxbl.^2 + Dybl.^2) ).*Dx(bl,:) + cosg*Dybl.*Dy(bl,:);
        L(bl,:) = (-ones(length(bl),1)-cosg*Dxbl./sqrt(1 + Dxbl.^2 + Dybl.^2)).*Dx(bl,:) + (-cosg*Dybl./sqrt(1 + Dxbl.^2 + Dybl.^2)).*Dy(bl,:);
        Nu(bl) = -Dxbl - cosg*sqrt(1 + Dxbl.^2 + Dybl.^2);

        Dxbr = Dxu(br);
        Dybr = Dyu(br);
        % L(br,:) = (cosg*Dxbr + sqrt(1 + Dxbr.^2 + Dybr.^2) ).*Dx(br,:) + cosg*Dybr.*Dy(br,:);
        % L(br,:) = (cosg*Dxbr - sqrt(1 + Dxbr.^2 + Dybr.^2) ).*Dx(br,:) + cosg*Dybr.*Dy(br,:);
        L(br,:) = (ones(length(br),1)-cosg*Dxbr./sqrt(1 + Dxbr.^2 + Dybr.^2)).*Dx(br,:) + (-cosg*Dybr./sqrt(1 + Dxbr.^2 + Dybr.^2)).*Dy(br,:);
        Nu(br) = Dxbr - cosg*sqrt(1 + Dxbr.^2 + Dybr.^2);

        Dxbur = Dxu(bur);
        Dybur = Dyu(bur);
        %L(bur,:) = (2*sqrt(2)*cosg*Dxbur - sqrt(1 + Dxbur.^2 + Dybur.^2) ).*Dx(bur,:) + (2*sqrt(2)*cosg*Dybur - sqrt(1 + Dxbur.^2 + Dybur.^2)).*Dy(bur,:);
        % L(bur,:) = (2*cosg*Dxbur - sqrt(1 + Dxbur.^2 + Dybur.^2) ).*Dx(bur,:) + (2*cosg*Dybur - sqrt(1 + Dxbur.^2 + Dybur.^2)).*Dy(bur,:);
        L(bur,:) = (1- 2*cosg*Dxbur./sqrt(1 + Dxbur.^2 + Dybur.^2)).*Dx(bur,:) + (1 -2*cosg*Dybur./sqrt(1 + Dxbur.^2 + Dybur.^2)).*Dy(bur,:);
        %Nu(bur) = Dxbur + Dybur - 2*sqrt(2)*cosg*sqrt(1 + Dxbur.^2 + Dybur.^2);
        Nu(bur) = Dxbur + Dybur - 2*cosg*sqrt(1 + Dxbur.^2 + Dybur.^2);

        Dxbul = Dxu(bul);
        Dybul = Dyu(bul);
        %L(bul,:) = (2*sqrt(2)*cosg*Dxbul + sqrt(1 + Dxbul.^2 + Dybul.^2) ).*Dx(bul,:) + (2*sqrt(2)*cosg*Dybul - sqrt(1 + Dxbul.^2 + Dybul.^2) ).*Dy(bul,:);
        % L(bul,:) = (2*cosg*Dxbul + sqrt(1 + Dxbul.^2 + Dybul.^2) ).*Dx(bul,:) + (2*cosg*Dybul - sqrt(1 + Dxbul.^2 + Dybul.^2) ).*Dy(bul,:);
        L(bul,:) = (-1- 2*cosg*Dxbul./sqrt(1 + Dxbul.^2 + Dybul.^2)).*Dx(bul,:) + (1 -2*cosg*Dybul./sqrt(1 + Dxbul.^2 + Dybul.^2)).*Dy(bul,:);
        %Nu(bul) = -Dxbul + Dybul - 2*sqrt(2)*cosg*sqrt(1 + Dxbul.^2 + Dybul.^2);
        Nu(bul) = -Dxbul + Dybul - 2*cosg*sqrt(1 + Dxbul.^2 + Dybul.^2);

        Dxbdr = Dxu(bdr);
        Dybdr = Dyu(bdr);
        %L(bdr,:) = (2*sqrt(2)*cosg*Dxbdr - sqrt(1 + Dxbdr.^2 + Dybdr.^2) ).*Dx(bdr,:) + (2*sqrt(2)*cosg*Dybdr + sqrt(1 + Dxbdr.^2 + Dybdr.^2)).*Dy(bdr,:);
        % L(bdr,:) = (2*cosg*Dxbdr - sqrt(1 + Dxbdr.^2 + Dybdr.^2) ).*Dx(bdr,:) + (2*cosg*Dybdr + sqrt(1 + Dxbdr.^2 + Dybdr.^2)).*Dy(bdr,:);
        L(bdr,:) = (1- 2*cosg*Dxbdr./sqrt(1 + Dxbdr.^2 + Dybdr.^2)).*Dx(bdr,:) + (-1 -2*cosg*Dybdr./sqrt(1 + Dxbdr.^2 + Dybdr.^2)).*Dy(bdr,:);
        %Nu(bdr) = Dxbdr - Dybdr - 2*sqrt(2)*cosg*sqrt(1 + Dxbdr.^2 + Dybdr.^2);
        Nu(bdr) = Dxbdr - Dybdr - 2*cosg*sqrt(1 + Dxbdr.^2 + Dybdr.^2);

        Dxbdl = Dxu(bdl);
        Dybdl = Dyu(bdl);
        %L(bdl,:) = (2*sqrt(2)*cosg*Dxbdl + sqrt(1 + Dxbdl.^2 + Dybdl.^2) ).*Dx(bdl,:) + (2*sqrt(2)*cosg*Dybdl + sqrt(1 + Dxbdl.^2 + Dybdl.^2)).*Dy(bdl,:);
        % L(bdl,:) = (2*cosg*Dxbdl + sqrt(1 + Dxbdl.^2 + Dybdl.^2) ).*Dx(bdl,:) + (2*cosg*Dybdl + sqrt(1 + Dxbdl.^2 + Dybdl.^2)).*Dy(bdl,:);
        L(bdl,:) = (-1- 2*cosg*Dxbdl./sqrt(1 + Dxbdl.^2 + Dybdl.^2)).*Dx(bdl,:) + (-1 -2*cosg*Dybdl./sqrt(1 + Dxbdl.^2 + Dybdl.^2)).*Dy(bdl,:);
        %Nu(bdl) = Dxbdl + Dybdl + 2*sqrt(2)*cosg*sqrt(1 + Dxbdl.^2 + Dybdl.^2);
        Nu(bdl) = -Dxbdl - Dybdl - 2*cosg*sqrt(1 + Dxbdl.^2 + Dybdl.^2);

        du = -L\Nu;
        new_res2 = (norm(du)) / (norm(u0)+ep);
        u0 = u0+du;
        count2 = count2 + 1;

        % uu0 = reshape(u0, N, N);
        % [xx,yy] = meshgrid(x,y);
        % 
        % uu0 = chebfun2(uu0);
        % figure(2)
        % plot(uu0)
        % pause
    end
    
    M = @(v) (Dxx*v) + (Dyy*v) + (Dxx*v).*(Dy*v).^2 - ...
              2*(Dx*v).*(Dy*v).*(Dxy*v) + (Dyy*v).*(Dx*v).^2 - ...
              ((1 + ((Dx*v).^2) + ((Dy*v).^2)).^(3/2).*kappa.*v);
    Mu = M(u0);
    Mui = Mu(inside);
    Nu(inside) = Mui;
    Dxu = Dx*u0;
    Dyu = Dy*u0;

    Dxbu = Dxu(bu);
    Dybu = Dyu(bu);
    Nu(bu) = Dybu - cosg*sqrt(1 + Dxbu.^2 + Dybu.^2);

    Dxbd = Dxu(bd);
    Dybd = Dyu(bd);
    Nu(bd) = -Dybd - cosg*sqrt(1 + Dxbd.^2 + Dybd.^2);

    Dxbl = Dxu(bl);
    Dybl = Dyu(bl);
    Nu(bl) = -Dxbl - cosg*sqrt(1 + Dxbl.^2 + Dybl.^2);

    Dxbr = Dxu(br);
    Dybr = Dyu(br);
    Nu(br) = Dxbr - cosg*sqrt(1 + Dxbr.^2 + Dybr.^2);

    Dxbur = Dxu(bur);
    Dybur = Dyu(bur);
    %Nu(bur) = Dxbur + Dybur - 2*sqrt(2)*cosg*sqrt(1 + Dxbur.^2 + Dybur.^2);
    Nu(bur) = Dxbur + Dybur - 2*cosg*sqrt(1 + Dxbur.^2 + Dybur.^2);
    
    Dxbul = Dxu(bul);
    Dybul = Dyu(bul);
    %Nu(bul) = -Dxbul + Dybul - 2*sqrt(2)*cosg*sqrt(1 + Dxbul.^2 + Dybul.^2);
    Nu(bul) = -Dxbul + Dybul - 2*cosg*sqrt(1 + Dxbul.^2 + Dybul.^2);
    
    Dxbdr = Dxu(bdr);
    Dybdr = Dyu(bdr);
    %Nu(bdr) = Dxbdr - Dybdr - 2*sqrt(2)*cosg*sqrt(1 + Dxbdr.^2 + Dybdr.^2);
    Nu(bdr) = Dxbdr - Dybdr - 2*cosg*sqrt(1 + Dxbdr.^2 + Dybdr.^2);
    
    Dxbdl = Dxu(bdl);
    Dybdl = Dyu(bdl);
    %Nu(bdl) = Dxbdl + Dybdl + 2*sqrt(2)*cosg*sqrt(1 + Dxbdl.^2 + Dybdl.^2);
    Nu(bdl) = -Dxbdl - Dybdl - 2*cosg*sqrt(1 + Dxbdl.^2 + Dybdl.^2);

    
    bvp_res = (norm(Nu)) / (norm(u0)+ep);
    
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

        inside = find(abs(xx)~=1 & abs(yy)~=1);
        % bu = find(abs(xx)<1 | yy==1);
        % bd = find(abs(xx)<1 | yy==-1);
        % bl = find((xx)==-1 | abs(yy)<1);
        % br = find((xx)==1 | abs(yy)<1);
        % bur = find(xx==1 | yy==1);
        % bul = find(xx==-1 | yy==1);
        % bdr = find(xx==1 | yy==-1);
        % bdl = find(xx==-1 | yy==-1);

        bu = find(abs(xx)<1 & yy==1);
        bd = find(abs(xx)<1 & yy==-1);
        bl = find((xx)==-1 & abs(yy)<1);
        br = find((xx)==1 & abs(yy)<1);
        bur = find(xx==1 & yy==1);
        bul = find(xx==-1 & yy==1);
        bdr = find(xx==1 & yy==-1);
        bdl = find(xx==-1 & yy==-1);

        M = @(v) (Dxx*v) + (Dyy*v) + (Dxx*v).*(Dy*v).^2 - ...
            2*(Dx*v).*(Dy*v).*(Dxy*v) + (Dyy*v).*(Dx*v).^2 - ...
            ((1 + ((Dx*v).^2) + ((Dy*v).^2)).^(3/2).*kappa.*v);
        Mu = M(u0);
        Mui = Mu(inside);
        Nu(inside) = Mui;
        
        Dxu = Dx*u0;
        Dyu = Dy*u0;

        Dxbu = Dxu(bu);
        Dybu = Dyu(bu);
        Nu(bu) = Dybu - cosg*sqrt(1 + Dxbu.^2 + Dybu.^2);

        Dxbd = Dxu(bd);
        Dybd = Dyu(bd);
        Nu(bd) = -Dybd - cosg*sqrt(1 + Dxbd.^2 + Dybd.^2);

        Dxbl = Dxu(bl);
        Dybl = Dyu(bl);
        Nu(bl) = -Dxbl - cosg*sqrt(1 + Dxbl.^2 + Dybl.^2);

        Dxbr = Dxu(br);
        Dybr = Dyu(br);
        Nu(br) = Dxbr - cosg*sqrt(1 + Dxbr.^2 + Dybr.^2);

        Dxbur = Dxu(bur);
        Dybur = Dyu(bur);
        %Nu(bur) = Dxbur + Dybur - 2*sqrt(2)*cosg*sqrt(1 + Dxbur.^2 + Dybur.^2);
        Nu(bur) = Dxbur + Dybur - 2*cosg*sqrt(1 + Dxbur.^2 + Dybur.^2);
        
        Dxbul = Dxu(bul);
        Dybul = Dyu(bul);
        %Nu(bul) = -Dxbul + Dybul - 2*sqrt(2)*cosg*sqrt(1 + Dxbul.^2 + Dybul.^2);
        Nu(bul) = -Dxbul + Dybul - 2*cosg*sqrt(1 + Dxbul.^2 + Dybul.^2);
        
        Dxbdr = Dxu(bdr);
        Dybdr = Dyu(bdr);
        %Nu(bdr) = Dxbdr - Dybdr - 2*sqrt(2)*cosg*sqrt(1 + Dxbdr.^2 + Dybdr.^2);
        Nu(bdr) = Dxbdr - Dybdr - 2*cosg*sqrt(1 + Dxbdr.^2 + Dybdr.^2);
        
        Dxbdl = Dxu(bdl);
        Dybdl = Dyu(bdl);
        %Nu(bdl) = Dxbdl + Dybdl + 2*sqrt(2)*cosg*sqrt(1 + Dxbdl.^2 + Dybdl.^2);
        Nu(bdl) = -Dxbdl - Dybdl - 2*cosg*sqrt(1 + Dxbdl.^2 + Dybdl.^2);
    end
    count1 = count1 + 1;

end

uu0 = reshape(u0, N, N);
[xx,yy] = meshgrid(x,y);

uu0 = chebfun2(uu0);
figure(3)
plot(uu0)
% axis equal
% figure(2)
% plot(uu0,'k-')