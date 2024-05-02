% Jonas Haug, Rachel Jewell, Ray Treinen, May 2024
%
% Compute minimal surfaces on the disk
% with Dirichlet data, defined to be function gamma.
%
% This function needs Chebfun installed to run: chebfun.org
% The dependencies on chebfun are restricted to the generation of the
% spectral differentiation matrices and plotting.

%% Physical parameters
g = @(t) 0.1*sin(2*t).^2;

%% Computational parameters
N = 50;
N2 = N/2;
M1 = 50;
M2 = M1/2;

new_tol = 1e-13;
bvp_tol = 1e-10;
ep = 1e-8;
MM = 100;

%% Computational building blocks
r = chebpts(N);
t = trigpts(M1, [-pi,pi]);

R = diag(r(N2+1:N));

Z = zeros(M2);
I = eye(M2);

[tt,rr] = meshgrid(t,r(N2+1:N));
rr = rr(:);
tt = tt(:);

D = diffmat(N, 1);
D2 = diffmat(N, 2);
D3 = D2(N2+1:N, N2:-1:1);
D4 = D2(N2+1:N,N2+1:N);
E3 = D(N2+1:N, N2:-1:1);
E4 = D(N2+1:N,N2+1:N);

Dr = kron([Z I; I Z],E3) + kron(eye(M1),E4);
Drr = kron([Z I; I Z],D3) + kron(eye(M1),D4);

D1t = diffmat(M1, 1, 'periodic', [-pi,pi]);
D2t = diffmat(M1, 2, 'periodic', [-pi,pi]);

Dth = kron(D1t,eye(N2));
Dthth = kron(D2t,eye(N2));
Drth = Dr * Dth;

% Initial guess
u0 = zeros(size(rr));

b = find(abs(rr)==1);
inside = find(abs(rr)~=1);

M = @(v) rr.*(Drr*v).*(rr.^2+(Dth*v).^2) + rr.*(Dthth*v).*(1+(Dr*v).^2) -...
    2*rr.*(Dr*v).*(Dth*v).*(Drth *v) + (Dr*v).*(rr.^2.*(1+(Dr*v).^2)+2*(Dth*v).^2);
Nu = zeros(size(u0));
Mu = M(u0);
Mui = Mu(inside);
Nu(inside) = Mui;
Nu(b) = u0(b) - g(tt(b));

%% Solving the problem
bvp_res = 1;
count1 = 0;
while((count1<MM) && (bvp_res > bvp_tol))
    new_res = 1;
    count2 = 0;
    z = spalloc(1, length(u0),1);
    while((count2 < MM) && (new_res > new_tol))

        M = @(v) rr.*(Drr*v).*(rr.^2+(Dth*v).^2) + rr.*(Dthth*v).*(1+(Dr*v).^2) -...
            2*rr.*(Dr*v).*(Dth*v).*(Drth *v) + (Dr*v).*(rr.^2.*(1+(Dr*v).^2)+2*(Dth*v).^2);

        F = @(v) rr.*((rr.^2) + (Dth*v).^2).*Drr + rr.*(1+(Dr*v).^2).*Dthth -...
            2*rr.*(Dr*v).*(Dth*v).*Drth + ...
            (2*rr.*(Dthth*v).*(Dr*v) - 2*rr.*(Dth*v).*(Drth*v) + rr.^2.*(1+3*(Dr*v).^2) + 2*(Dth*v).^2).*Dr +...
            (4*(Dr*v).*(Dth*v) - 2*rr.*(Dr*v).*(Drth*v) + 2*rr.*(Dth*v).*(Drr*v)).*Dth;

        Mu = M(u0);
        Mui = Mu(inside);
        Nu(inside) = Mui;
        Nu(b) = u0(b)-g(tt(b));
        
        L = F(u0);

        for ii=1:length(b)
            L(b(ii),:) = z;
            L(b(ii),b(ii)) = 1;
        end

        du = -L\Nu;
        
        new_res = norm(du)/(norm(u0) + ep);
        u0 = u0 + du;
        count2 = count2+1;

    end

    M = @(v) rr.*(Drr*v).*(rr.^2+(Dth*v).^2) + rr.*(Dthth*v).*(1+(Dr*v).^2) -...
        2*rr.*(Dr*v).*(Dth*v).*(Drth *v) + (Dr*v).*(rr.^2.*(1+(Dr*v).^2)+2*(Dth*v).^2);
    Mu = M(u0);
    Mui = Mu(inside);
    Nu(inside) = Mui;
    Nu(b) = u0(b) - g(tt(b));

    bvp_res = norm(Nu)/(norm(u0) + ep);

    if (bvp_res > bvp_tol)
        %         uu0 = reshape(u0,N,N);
        %         uu0 = chebfun2(uu0);

        M1 = M1+10;
        N = N+10;
        N2 = N/2;
        M2 = M1/2;

        r = chebpts(N);
        t = trigpts(M1, [0;2*pi]);

        R = diag(r(N2+1:N));

        Z = zeros(M2);
        I = eye(M2);

        [tt,rr] = meshgrid(t,r(N2+1:N));
        rr = rr(:);
        tt = tt(:);

        D = diffmat(N, 1);
        D2 = diffmat(N, 2);
        D3 = D2(N2+1:N, N2:-1:1);
        D4 = D2(N2+1:N,N2+1:N);
        E3 = D(N2+1:N, N2:-1:1);
        E4 = D(N2+1:N,N2+1:N);

        Dr = kron([Z I; I Z],E3) + kron(eye(M1),E4);
        Drr = kron([Z I; I Z],D3) + kron(eye(M1),D4);

        D1t = diffmat(M1, 1, 'periodic', [-pi,pi]);
        D2t = diffmat(M1, 2, 'periodic', [-pi,pi]);

        Dth = kron(D1t,eye(N2));
        Dthth = kron(D2t,eye(N2));
        Drth = Dr * Dth;

        u0 = zeros(size(rr));
        b = find(abs(rr)==1);
        inside = find(abs(rr)~=1);

        M = @(v) rr.*(Drr*v).*(rr.^2+(Dth*v).^2) + rr.*(Dthth*v).*(1+(Dr*v).^2) -...
            2*rr.*(Dr*v).*(Dth*v).*(Drth *v) + (Dr*v).*(rr.^2.*(1+(Dr*v).^2)+2*(Dth*v).^2);
        Nu = zeros(size(u0));
        Mu = M(u0);
        Mui = Mu(inside);
        Nu(inside) = Mui;
        Nu(b) = u0(b) - g(tt(b));
    end
    count1 = count1 + 1;
end

%% Plotting

uu0 = reshape(u0,N2,M1);
uY = uu0;
uu0 = uu0(:,[M1 1:M1]);
[tt,rr] = meshgrid(t([M1 1:M1]),r(N2+1:N));
[xx,yy] = pol2cart(tt,rr);

Y = diskfun(uY);

figure(2)
surf(Y)
xlabel('X', 'FontWeight', 'bold')
ylabel('Y', 'FontWeight', 'bold')
zlabel('U', 'FontWeight', 'bold')
fontsize("increase")
% axis equal

figure(3)
contour(Y)
xlabel('X', 'FontWeight', 'bold')
ylabel('Y', 'FontWeight', 'bold')
fontsize("increase")
axis equal

