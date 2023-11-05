a = 1; 
b = 2;

N = 50;
M1 = 80;

new_tol = 1e-13;
bvp_tol = 1e-10;
ep = 1e-8;
MM = 100;

r = chebpts(N, [a;b]);
t = trigpts(M1, [-pi,pi]);
R = diag(r);

h = 1;
ha = h;
hb = 0;

[tt,rr] = meshgrid(t,r);
rr = rr(:);
tt = tt(:);

D = diffmat(N, 1, [a;b]);
D2 = diffmat(N, 2, [a;b]);

Dr = kron(eye(M1),D);
Drr = kron(eye(M1),D2);

D1t = diffmat(M1, 1, 'periodic', [-pi,pi]);
D2t = diffmat(M1, 2, 'periodic', [-pi,pi]);

Dth = kron(D1t,eye(N));
Dthth = kron(D2t,eye(N));
Drth = Dr * Dth;

u0 = zeros(size(rr));
ba = find(rr==a);
bb = find(rr==b);
inside = find((rr~=a)&(rr~=b));
g = @(t) 0.1*sin(2*t).^2;

M = @(v) rr.*(Drr*v).*(rr.^2 + (Dth*v).^2) + rr.*(Dthth*v).*(1 + (Dr*v).^2) - 2*rr.*(Dr*v).*(Dth*v).*(Drth *v) + (Dr*v).*(rr.^2.*(1 + (Dr*v).^2) + 2*(Dth*v).^2);
Nu = zeros(size(u0));
Mu = M(u0);
Mui = Mu(inside);
Nu(inside) = Mui;

Nu(ba) = u0(ba) - ha;
Nu(bb) = u0(bb) - hb;

bvp_res = 1;
count1 = 0;
while((count1<MM) && (bvp_res > bvp_tol))
    new_res = 1;
    count2 = 0;
    z = spalloc(1, length(u0),1);
    while((count2 < MM) && (new_res > new_tol))

        M = @(v) rr.*(Drr*v).*(rr.^2 + (Dth*v).^2) + rr.*(Dthth*v).*(1 + (Dr*v).^2) - 2*rr.*(Dr*v).*(Dth*v).*(Drth *v) + (Dr*v).*(rr.^2.*(1 + (Dr*v).^2) + 2*(Dth*v).^2);

        F = @(v) rr.*((rr.^2) + (Dth*v).^2).*Drr + rr.*(1 + (Dr*v).^2).*Dthth - 2*rr.*(Dr*v).*(Dth*v).*Drth + (2*rr.*(Dthth*v).*(Dr*v) - 2*rr.*(Dth*v).*(Drth*v) + rr.^2.*(1+3*(Dr*v).^2) + 2*(Dth*v).^2).*Dr + (4*(Dr*v).*(Dth*v) - 2*rr.*(Dr*v).*(Drth*v)  + 2*rr.*(Dth*v).*(Drr*v)).*Dth;

        Mu = M(u0);
        Mui = Mu(inside);
        Nu(inside) = Mui;
        Nu(ba) = u0(ba) - ha;
        Nu(bb) = u0(bb) - hb;

        L = F(u0);

        for ii=1:length(ba)
            L(ba(ii),:) = z;
            L(ba(ii),ba(ii)) = 1;
        end

        for ii=1:length(bb)
            L(bb(ii),:) = z;
            L(bb(ii),bb(ii)) = 1;
        end

        du = -L\Nu;
        new_res = norm(du)/(norm(u0)+ep);
        u0 = u0+du;
        count2 = count2+1;

    end

    M = @(v) rr.*(Drr*v).*(rr.^2 + (Dth*v).^2) + rr.*(Dthth*v).*(1 + (Dr*v).^2) - 2*rr.*(Dr*v).*(Dth*v).*(Drth *v) + (Dr*v).*(rr.^2.*(1 + (Dr*v).^2) + 2*(Dth*v).^2);
    Mu = M(u0);
    Mui = Mu(inside);
    Nu(inside) = Mui;
    Nu(ba) = u0(ba) - ha;
    Nu(bb) = u0(bb) - hb;

    bvp_res = norm(Nu)/(norm(u0)+ep);

    if (bvp_res > bvp_tol)
        %         uu0 = reshape(u0,N,N);
        %         uu0 = chebfun2(uu0);

        M1 = M1+10;
        N = N+10;

        r = chebpts(N);
        t = trigpts(M1, [0;2*pi]);

        R = diag(r);

        [tt,rr] = meshgrid(t,r);
        rr = rr(:);
        tt = tt(:);

        D = diffmat(N, 1, [a;b]);
        D2 = diffmat(N, 2, [a;b]);

        Dr = kron(eye(M1),D);
        Drr = kron(eye(M1),D2);

        D1t = diffmat(M1, 1, 'periodic', [-pi,pi]);
        D2t = diffmat(M1, 2, 'periodic', [-pi,pi]);

        Dth = kron(D1t,eye(N));
        Dthth = kron(D2t,eye(N));
        Drth = Dr * Dth;

        % R = diag(r);

        u0 = zeros(size(rr));
        ba = find(rr==a);
        bb = find(rr==b);
        inside = find((rr~=a)&(rr~=b));
        % g = @(t) 0.1*sin(2*t).^2;
        %g = @(t) 0;

        M = @(v) rr.*(Drr*v).*(rr.^2 + (Dth*v).^2) + rr.*(Dthth*v).*(1 + (Dr*v).^2) - 2*rr.*(Dr*v).*(Dth*v).*(Drth *v) + (Dr*v).*(rr.^2.*(1 + (Dr*v).^2) + 2*(Dth*v).^2);
        Nu = zeros(size(u0));
        Mu = M(u0);
        Mui = Mu(inside);
        Nu(inside) = Mui;
        Nu(ba) = u0(ba) - ha;
        Nu(bb) = u0(bb) - hb;
    end
    count1 = count1 + 1;
end
% length(u0) - N*M1/2
uu0 = reshape(u0,N,M1);
% isreal(uu0)
% uu0 = real(uu0);
uY = uu0;
uu0 = uu0(:,[M1 1:M1]);
[tt,rr] = meshgrid(t([M1 1:M1]),r);
[xx,yy] = pol2cart(tt,rr);

surf(xx,yy,uu0)
xlabel x
ylabel y
zlabel u
