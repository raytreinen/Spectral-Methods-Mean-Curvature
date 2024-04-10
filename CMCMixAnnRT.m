% Jonas Haug, Rachel Jewell, Ray Treinen, March 2024
% 
% Compute constant mean curvature surfaces on the annulus with mixed
% Dirichlet and Neumann data. 
% For the inner radius, r = a and the data is defined to be function 
% ha or agamma for Dirichlet or Neumann, respectively.
% For the outter radius, r = b and the data is defined to be function
% hb or bgamma for Dirichlet or Neumann, respectively.
% Kappa is set according to the physical problem.
%
% This function needs Chebfun installed to run: chebfun.org
% The dependencies on chebfun are restricted to the generation of the
% spectral differentiation matrices and plotting.

%% Physical parameters
a = 1; 
b = 2;

lambda = -0.85;

% Neumann boundary functions
agamma = @(t) pi/2 + 0.2*sin(4*t);
% bgamma = @(t) pi/2 + 0.1*sin(4*t);

% Dirichlet boundary functions
% ha = @(t) 0.5 + 0.1*sin(2*t).^2;
hb = @(t) 0.5 - 0.1*cos(2*t).^2;

%% Computational parameters
N = 50;
M1 = 80;

new_tol = 1e-13;
bvp_tol = 1e-10;
ep = 1e-8;
MM = 100;

%% Computational building blocks
r = chebpts(N, [a;b]);
t = trigpts(M1, [-pi,pi]);
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

% Neumann boundary functions
acg = @(t) cos(agamma(t));
acgc = chebfun(acg,[0,2*pi]);
aic = sum(acgc,0,2*pi);
alambda = aic/pi;

% bcg = @(t) cos(bgamma(t));
% bcgc = chebfun(bcg,[0,2*pi]);
% bic = sum(bcgc,0,2*pi);
% blambda = bic/pi;

% Initial Guess
u0 = ones(size(rr));
ba = find(rr==a);
bb = find(rr==b);
inside = find((rr~=a)&(rr~=b));

M = @(v) rr.*(Drr*v).*(rr.^2+(Dth*v).^2) + rr.*(Dthth*v).*(1+(Dr*v).^2) -...
    2*rr.*(Dr*v).*(Dth*v).*(Drth *v) + (Dr*v).*(rr.^2.*(1+(Dr*v).^2)+2*(Dth*v).^2) -...
    lambda*(rr.^2.*(1 + (Dr*v).^2) + (Dth*v).^2).^(3/2);Nu = zeros(size(u0));
Mu = M(u0);
Mui = Mu(inside);
Nu(inside) = Mui;

% For Neumann data
Drba = Dr*u0;
Drba = Drba(ba);
Dthba = Dth*u0;
Dthba = Dthba(ba);
Nu(ba) = rr(ba).*Drba + cos(agamma(tt(ba))).*sqrt(rr(ba).^2.*(1 + Drba.^2) + Dthba.^2);

% Drbb = Dr*u0;
% Drbb = Drbb(bb);
% Dthbb = Dth*u0;
% Dthbb = Dthbb(bb);
% Nu(bb) = rr(bb).*Drbb - cos(bgamma(tt(bb))).*sqrt(rr(bb).^2.*(1 + Drbb.^2) + Dthbb.^2);

% For Dirichlet data
% Nu(ba) = u0(ba) - ha(tt(ba));
Nu(bb) = u0(bb) - hb(tt(bb));

%% Solving the problem
bvp_res = 1;
count1 = 0;
while((count1<MM) && (bvp_res > bvp_tol))
    new_res = 1;
    count2 = 0;
    z = spalloc(1, length(u0),1);
    while((count2 < MM) && (new_res > new_tol))

        M = @(v) rr.*(Drr*v).*(rr.^2+(Dth*v).^2) + rr.*(Dthth*v).*(1+(Dr*v).^2) -...
            2*rr.*(Dr*v).*(Dth*v).*(Drth *v) + (Dr*v).*(rr.^2.*(1+(Dr*v).^2)+2*(Dth*v).^2) -...
            lambda*(rr.^2.*(1 + (Dr*v).^2) + (Dth*v).^2).^(3/2);

        F = @(v) rr.*((rr.^2) + (Dth*v).^2).*Drr + rr.*(1 + (Dr*v).^2).*Dthth -...
            2*rr.*(Dr*v).*(Dth*v).*Drth + ...
            (2*rr.*(Dthth*v).*(Dr*v) - 2*rr.*(Dth*v).*(Drth*v) + rr.^2.*(1 + 3*(Dr*v).^2) + 2*(Dth*v).^2).*Dr +...
            (4*(Dr*v).*(Dth*v) - 2*rr.*(Dr*v).*(Drth*v) + 2*rr.*(Dth*v).*(Drr*v)).*Dth -...
            3*lambda*sqrt((rr.^2).*(1 + (Dr*v)) + (Dth.^2)).*(rr.^2.*(Dr*v).*Dr + (Dth*v).*Dth);

        Mu = M(u0);
        Mui = Mu(inside);
        Nu(inside) = Mui;
        
        % Neumann boundary
        Drba = Dr*u0;
        Drba = Drba(ba);
        Dthba = Dth*u0;
        Dthba = Dthba(ba);
        Nu(ba) = rr(ba).*Drba + cos(agamma(tt(ba))).*sqrt(rr(ba).^2.*(1 + Drba.^2) + Dthba.^2);

        % Drbb = Dr*u0;
        % Drbb = Drbb(bb);
        % Dthbb = Dth*u0;
        % Dthbb = Dthbb(bb);
        % Nu(bb) = rr(bb).*Drbb - cos(bgamma(tt(bb))).*sqrt(rr(bb).^2.*(1 + Drbb.^2) + Dthbb.^2);

        % Dirichlet boundary
        % Nu(ba) = u0(ba) - ha(tt(ba));
        Nu(bb) = u0(bb) - hb(tt(bb));

        L = F(u0);

        % Neumann Frechet
        L(ba,:) = (rr(ba) + cos(agamma(tt(ba))).*rr(ba).^2.*Drba./sqrt(rr(ba).^2.*(1 + Drba.^2) + Dthba.^2) ).*Dr(ba,:) + ...
            (cos(agamma(tt(ba))).*Dthba./sqrt(rr(ba).^2.*(1 + Drba.^2) + Dthba.^2)).*Dth(ba,:);
        % L(bb,:) = (rr(bb) - cos(bgamma(tt(bb))).*rr(bb).^2.*Drbb./sqrt(rr(bb).^2.*(1 + Drbb.^2) + Dthbb.^2) ).*Dr(bb,:) - (cos(bgamma(tt(bb))).*Dthbb./sqrt(rr(bb).^2.*(1 + Drbb.^2) + Dthbb.^2)).*Dth(bb,:);

        % Dirichlet Frechet
        % for ii=1:length(ba)
        %     L(ba(ii),:) = z;
        %     L(ba(ii),ba(ii)) = 1;
        % end

        for ii=1:length(bb)
            L(bb(ii),:) = z;
            L(bb(ii),bb(ii)) = 1;
        end

        du = -L\Nu;

        new_res = norm(du)/(norm(u0)+ep);
        u0 = u0 + du;
        
        % uY = reshape(u0,N2,M1);
        % Y = diskfun(uY);
        % figure(20)
        % surf(Y)
        % pause
        
        count2 = count2 + 1;

    end

    M = @(v) rr.*(Drr*v).*(rr.^2+(Dth*v).^2) + rr.*(Dthth*v).*(1+(Dr*v).^2) -...
        2*rr.*(Dr*v).*(Dth*v).*(Drth *v) + (Dr*v).*(rr.^2.*(1+(Dr*v).^2)+2*(Dth*v).^2) -...
        lambda*(rr.^2.*(1 + (Dr*v).^2) + (Dth*v).^2).^(3/2);
    Mu = M(u0);
    Mui = Mu(inside);
    Nu(inside) = Mui;
    
    % Neuman
    Drba = Dr*u0;
    Drba = Drba(ba);
    Dthba = Dth*u0;
    Dthba = Dthba(ba);
    Nu(ba) = rr(ba).*Drba + cos(agamma(tt(ba))).*sqrt(rr(ba).^2.*(1 + Drba.^2) + Dthba.^2);

    % Drbb = Dr*u0;
    % Drbb = Drbb(bb);
    % Dthbb = Dth*u0;
    % Dthbb = Dthbb(bb);
    % Nu(bb) = rr(bb).*Drbb - cos(bgamma(tt(bb))).*sqrt(rr(bb).^2.*(1 + Drbb.^2) + Dthbb.^2);

    % Dirichlet
    % Nu(ba) = u0(ba) - ha(tt(ba));    
    Nu(bb) = u0(bb) - hb(tt(bb));

    bvp_res = norm(Nu)/(norm(u0)+ep);

    if (bvp_res > bvp_tol)
        %         uu0 = reshape(u0,N,N);
        %         uu0 = chebfun2(uu0);

        M1 = M1+10;
        N = N+10;

        r = chebpts(N, [a;b]);
        t = trigpts(M1, [-pi,pi]);
        R = diag(r);

        [tt,rr] = meshgrid(t,r);
        rr = rr(:);
        tt = tt(:);

        D = diffmat(N, 1);
        D2 = diffmat(N, 2);

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

        M = @(v) rr.*(Drr*v).*(rr.^2+(Dth*v).^2) + rr.*(Dthth*v).*(1+(Dr*v).^2) -...
            2*rr.*(Dr*v).*(Dth*v).*(Drth *v) + (Dr*v).*(rr.^2.*(1+(Dr*v).^2)+2*(Dth*v).^2) -...
            lambda*(rr.^2.*(1 + (Dr*v).^2) + (Dth*v).^2).^(3/2);
        Nu = zeros(size(u0));
        Mu = M(u0);
        Mui = Mu(inside);
        Nu(inside) = Mui;
        
        Drba = Dr*u0;
        Drba = Drba(ba);
        Dthba = Dth*u0;
        Dthba = Dthba(ba);
        Nu(ba) = rr(ba).*Drba + cos(agamma(tt(ba))).*sqrt(rr(ba).^2.*(1 + Drba.^2) + Dthba.^2);

        % Drbb = Dr*u0;
        % Drbb = Drbb(bb);
        % Dthbb = Dth*u0;
        % Dthbb = Dthbb(bb);
        % Nu(bb) = rr(bb).*Drbb - cos(bgamma(tt(bb))).*sqrt(rr(bb).^2.*(1 + Drbb.^2) + Dthbb.^2);

        % Nu(ba) = u0(ba) - ha(tt(ba));
        Nu(bb) = u0(bb) - hb(tt(bb));
    end
    count1 = count1 + 1;
end

%% Plotting
% length(u0) - N*M1/2
uu0 = reshape(u0,N,M1);
uY = uu0;
uu0 = uu0(:,[M1 1:M1]);
[tt,rr] = meshgrid(t([M1 1:M1]),r);
[xx,yy] = pol2cart(tt,rr);

figure(1)
surf(xx,yy,uu0)
xlabel('X', 'FontWeight', 'bold')
ylabel('Y', 'FontWeight', 'bold')
zlabel('U', 'FontWeight', 'bold')
fontsize("increase")
% axis equal

figure(3)
contour(xx,yy,uu0)
xlabel('X', 'FontWeight', 'bold')
ylabel('Y', 'FontWeight', 'bold')
fontsize("increase")
axis equal
hold on
plot(a*cos(t([M1 1:M1])), a*sin(t([M1 1:M1])), 'k')
plot(b*cos(t([M1 1:M1])), b*sin(t([M1 1:M1])), 'k')
