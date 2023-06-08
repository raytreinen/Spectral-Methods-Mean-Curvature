function [] = symmetric_capillary_disc_3z()
% Ray Treinen, April 2022
%
% Compute capillary surfaces
% with inclination angles psia at radius r = a and psib at radius r = b
%
% The output can be set to [R, U, Psi, ell] to return the solution to the
% problem.  The default is set to merely plot the solution.
%
% This function needs Chebfun installed to run: chebfun.org
% The dependencies on chebfun are restricted to the generation of the
% spectral differentiation matrices and plotting.
%
% If this file is used as a function, replace the first line with
% function [R, U, Psi ,ell, n, res_bvp] = symmetric_capillary_disc_damped()
% and delete (or comment) the physical parameters block as well as the figure plotting
% at the end of the program.

%% physical parameters
% this block can be commented out if the program is used inside of a larger
% program and these values are passed as inputs to this program.
kappa = 1;
b = 40;
% d is only used if b > 1
d = 5;
psib_actual = 8*pi/8;

%% Computational parameters
% computational domain
X = [-1,1];

% maximum loop counts and the tolerances for each loop
max_iter_newton = 100000;
max_iter_bvp = 10000;
tol_newton = 1e-13;
tol_bvp = 1e-12;

% initialize the number of points.
% this is chosen so that the singularity in the ODE is avoided.
k = 7;
n = 2*k + 1;
n1 = k;
n2 = 2*k + 1;
n3 = k;

%% Initial guesses

% assign to psia and psib the actual values if they are less than pi/2 in
% magnitude, and if not, then replace them with +- pi/2
psib = sign(psib_actual)*min(abs(psib_actual), pi/2);

if sqrt(kappa)*b >= 1
    % compute the Chebyshev points (called tau in the paper)
    s1 = chebpts(n1);
    s2 = chebpts(n2);
    s3 = chebpts(n3);

    % compute initial guesses for psi, R and U.  These may be replaced below.
    Psi01 = @(s1) (1 - s1)*(-psib)/2;
    Psi02 = zeros(n2,1);
    Psi03 = @(s3) (1 + s3)*psib/2;

    u02 = zeros(n2,1);
    R2 = (b-d)*s2;
    ell02 = 2*(b-d);

    if (psib >= 0)
        R1 = ((d-b) - (-b))/(sin(-(-psib)) + sin(0));
        R3 = (b - (b-d))/(sin(0) + sin(psib));
        ell01 = R1*(-(-psib));
        ell03 = R3*(psib);

        R01 = ((-b)+(d-b))/2 + R1*sin(Psi01(s1));
        R02 = R2;
        R03 = ((b-d)+b)/2 + R3*sin(Psi03(s3));
        U01 = 1/R1 + R1*(1-cos(Psi01(s1)));
        U02 = u02;
        U03 = 1/R3  + R3*(1-cos(Psi03(s3)));
        P01 = Psi01(s1);
        P02 = Psi02;
        P03 = Psi03(s3);
        ell01;
        ell02;
        ell03;

        v = [R01; R02; R03; U01; U02; U03; P01; P02; P03; ell01; ell02; ell03];

    else
        R1 = ((d-b) + b)/(sin(-psib) + sin(0));
        R3 = (b - (b-d))/(sin(0) + sin(-psib));
        ell01 = R1*(-psib);
        ell03 = R3*(-psib);

        R01 = ((-b)+(d-b))/2 + R1*sin(-Psi01(s1));
        R02 = R2;
        R03 = ((b-d)+b)/2 + R3*sin(-Psi03(s3));
        U01 = -1/R1 - R1*(1-cos(Psi01(s1)));
        U02 = u02;
        U03 = -1/R3  - R3*(1-cos(Psi03(s3)));
        P01 = Psi01(s1);
        P02 = Psi02;
        P03 = Psi03(s3);
        ell01;
        ell02;
        ell03;

        v = [R01; R02; R03; U01; U02; U03; P01; P02; P03; ell01; ell02; ell03];

    end
    %% Solving the problem
    % the cases depend on if the solution is a graph over a base domain

    if (abs(psib_actual) > pi/2)
        [v, n1, n2, n3] = cheb_engine_large(v, n1, n2, n3);

        psi_vec = linspace(sign(psib_actual)*pi/2,psib_actual,11)';
        psi_vec = psi_vec(2:end);

        for i = 1:10
            psib = psi_vec(i);
            [v, n1, n2, n3] = cheb_engine_large(v, n1, n2, n3);
        end
    else
        [v, n1, n2, n3] = cheb_engine_large(v, n1, n2, n3);
    end

else
    % compute the Chebyshev points (called tau in the paper)
    ss = chebpts(n);

    % compute initial guesses for psi, R and U.  These may be replaced below.
    Psi0 = @(s) s*psib;

    gma = pi/2 - abs(psib);
    % u0 = 2*cos(gma)/(kappa*b) - b/cos(gma) + (2*b/3)*(1-sin(gma)^3)/cos(gma)^3;
    u0 = 2*cos(gma)/(kappa*b) - b*cos(gma)*(1 + 2*sin(gma))/(3*(1 + sin(gma))^2);
    % average of the two radii from the upper and lower estimates
    R = 0.5*(2/(u0*kappa) + abs(b/sin(psib)));

    ell0 = abs(R*psib);
    v = [sign(psib)*R*sin(Psi0(ss));sign(psib)*( R*(1 - cos(Psi0(ss))) + u0 );Psi0(ss); ell0];

    %% Solving the problem
    % the cases depend on if the solution is a graph over a base domain
    if (abs(psib_actual) > pi/2)
        [v, ~] = cheb_engine_small(v, n);

        n = (length(v) - 1)/3;
        psi_vec = linspace(sign(psib_actual)*pi/2,psib_actual,11)';
        psi_vec = psi_vec(2:end);

        for i = 1:10
            psib = psi_vec(i);
            [v, n] = cheb_engine_small(v, n);
        end
    else
        [v, n] = cheb_engine_small(v, n);
    end
end

%% the main computational engine of the file
    function [v, n1, n2, n3] = cheb_engine_large(v, n1, n2, n3)

        % intialize the residual
        res_bvp1 = 1;
        res_bvp2 = 1;
        res_bvp3 = 1;

        while (res_bvp1 > tol_bvp)||(res_bvp2 > tol_bvp)||(res_bvp3 > tol_bvp)

            % initialize the differential operator components
            %
            % D0 and D1 are spectral operators corresponding to a
            % downsampling matrix and a differentiation matrix,
            % respectively
            D0 = diffmat([n1-1 n1],0,X);
            D1 = diffmat([n1-1 n1],1,X);
            Z0 = sparse(n1-1,n1);
            D00 = diffmat([n2-1 n2],0,X);
            D10 = diffmat([n2-1 n2],1,X);
            Z00 = sparse(n2-1,n2);
            D000 = diffmat([n3-1 n3],0,X);
            D100 = diffmat([n3-1 n3],1,X);
            Z000 = sparse(n3-1,n3);
            D01 = spalloc(n1-1, 3*(n1 + n2 + n3 + 1), n1*(n1-1));
            D02 = spalloc(n2-1, 3*(n1 + n2 + n3 + 1), n2*(n2-1));
            D03 = spalloc(n3-1, 3*(n1 + n2 + n3 + 1), n3*(n3-1));
            D04 = D01;
            D05 = D02;
            D06 = D03;
            D07 = D01;
            D08 = D02;
            D09 = D03;
            D11 = D01;
            D12 = D02;
            D13 = D03;
            D14 = D01;
            D15 = D02;
            D16 = D03;
            D17 = D01;
            D18 = D02;
            D19 = D03;
            D01(1:n1-1, 1:n1) = D0;
            D11(1:n1-1, 1:n1) = D1;
            D02(1:n2-1, n1+1:n1+n2) = D00;
            D12(1:n2-1, n1+1:n1+n2) = D10;
            D03(1:n3-1, n1+n2+1:n1+n2+n3) = D000;
            D13(1:n3-1, n1+n2+1:n1+n2+n3) = D100;
            D04(1:n1-1, n1+n2+n3+1:2*n1+n2+n3) = D0;
            D14(1:n1-1, n1+n2+n3+1:2*n1+n2+n3) = D1;
            D05(1:n2-1, 2*n1+n2+n3+1:2*n1+2*n2+n3) = D00;
            D15(1:n2-1, 2*n1+n2+n3+1:2*n1+2*n2+n3) = D10;
            D06(1:n3-1, 2*n1+2*n2+n3+1:2*n1+2*n2+2*n3) = D000;
            D16(1:n3-1, 2*n1+2*n2+n3+1:2*n1+2*n2+2*n3) = D100;
            D07(1:n1-1, 2*n1+2*n2+2*n3+1:3*n1+2*n2+2*n3) = D0;
            D17(1:n1-1, 2*n1+2*n2+2*n3+1:3*n1+2*n2+2*n3) = D1;
            D08(1:n2-1, 3*n1+2*n2+2*n3+1:3*n1+3*n2+2*n3) = D00;
            D18(1:n2-1, 3*n1+2*n2+2*n3+1:3*n1+3*n2+2*n3) = D10;
            D09(1:n3-1, 3*n1+3*n2+2*n3+1:3*n1+3*n2+3*n3) = D000;
            D19(1:n3-1, 3*n1+3*n2+2*n3+1:3*n1+3*n2+3*n3) = D100;

            DN1_01 = spalloc(n1-1, 3*n1+1, n1*(n1-1));
            DN2_02 = spalloc(n2-1, 3*n2+1, n2*(n2-1));
            DN3_03 = spalloc(n3-1, 3*n3+1, n3*(n3-1));
            DN1_04 = DN1_01;
            DN2_05 = DN2_02;
            DN3_06 = DN3_03;
            DN1_07 = DN1_01;
            DN2_08 = DN2_02;
            DN3_09 = DN3_03;
            DN1_11 = DN1_01;
            DN2_12 = DN2_02;
            DN3_13 = DN3_03;
            DN1_14 = DN1_01;
            DN2_15 = DN2_02;
            DN3_16 = DN3_03;
            DN1_17 = DN1_01;
            DN2_18 = DN2_02;
            DN3_19 = DN3_03;
            DN1_01(1:n1-1, 1:n1) = D0;
            DN1_11(1:n1-1, 1:n1) = D1;
            DN2_02(1:n2-1, 1:n2) = D00;
            DN2_12(1:n2-1, 1:n2) = D10;
            DN3_03(1:n3-1, 1:n3) = D000;
            DN3_13(1:n3-1, 1:n3) = D100;
            DN1_04(1:n1-1, n1+1:2*n1) = D0;
            DN1_14(1:n1-1, n1+1:2*n1) = D1;
            DN2_05(1:n2-1, n2+1:2*n2) = D00;
            DN2_15(1:n2-1, n2+1:2*n2) = D10;
            DN3_06(1:n3-1, n3+1:2*n3) = D000;
            DN3_16(1:n3-1, n3+1:2*n3) = D100;
            DN1_07(1:n1-1, 2*n1+1:3*n1) = D0;
            DN1_17(1:n1-1, 2*n1+1:3*n1) = D1;
            DN2_08(1:n2-1, 2*n2+1:3*n2) = D00;
            DN2_18(1:n2-1, 2*n2+1:3*n2) = D10;
            DN3_09(1:n3-1, 2*n3+1:3*n3) = D000;
            DN3_19(1:n3-1, 2*n3+1:3*n3) = D100;


            % Evaluate the computational vector to check the boundary
            % conditions
            dT1n1 = sparse(1,3*(n1 + n2 + n3 + 1));
            dT1p1 = dT1n1;
            dT2n1 = dT1n1;
            dT2p1 = dT1n1;
            dT3n1 = dT1n1;
            dT3p1 = dT1n1;
            dT4n1 = dT1n1;
            dT4p1 = dT1n1;
            dT5n1 = dT1n1;
            dT5p1 = dT1n1;
            dT6n1 = dT1n1;
            dT6p1 = dT1n1;
            dT7n1 = dT1n1;
            dT7p1 = dT1n1;
            dT8n1 = dT1n1;
            dT8p1 = dT1n1;
            dT9n1 = dT1n1;
            dT9p1 = dT1n1;
            dT1n1(1) = 1;
            dT1p1(n1) = 1;
            dT2n1(n1+1) = 1;
            dT2p1(n1+n2) = 1;
            dT3n1(n1+n2+1) = 1;
            dT3p1(n1+n2+n3) = 1;
            dT4n1(n1+n2+n3+1) = 1;
            dT4p1(2*n1+n2+n3) = 1;
            dT5n1(2*n1+n2+n3+1) = 1;
            dT5p1(2*n1+2*n2+n3) = 1;
            dT6n1(2*n1+2*n2+n3+1) = 1;
            dT6p1(2*n1+2*n2+2*n3) = 1;
            dT7n1(2*n1+2*n2+2*n3+1) = 1;
            dT7p1(3*n1+2*n2+2*n3) = 1;
            dT8n1(3*n1+2*n2+2*n3+1) = 1;
            dT8p1(3*n1+3*n2+2*n3) = 1;
            dT9n1(3*n1+3*n2+2*n3+1) = 1;
            dT9p1(end-3) = 1;

            dTN1_1n1 = sparse(1,3*n1+1);
            dTN1_1p1 = dTN1_1n1;
            dTN2_2n1 = sparse(1,3*n2+1);
            dTN2_2p1 = dTN2_2n1;
            dTN3_3n1 = sparse(1,3*n3+1);
            dTN3_3p1 = dTN3_3n1;
            dTN1_4n1 = dTN1_1n1;
            dTN1_4p1 = dTN1_1n1;
            dTN2_5n1 = dTN2_2n1;
            dTN2_5p1 = dTN2_2n1;
            dTN3_6n1 = dTN3_3n1;
            dTN3_6p1 = dTN3_3n1;
            dTN1_7n1 = dTN1_1n1;
            dTN1_7p1 = dTN1_1n1;
            dTN2_8n1 = dTN2_2n1;
            dTN2_8p1 = dTN2_2n1;
            dTN3_9n1 = dTN3_3n1;
            dTN3_9p1 = dTN3_3n1;
            dTN1_1n1(1) = 1;
            dTN1_1p1(n1) = 1;
            dTN2_2n1(1) = 1;
            dTN2_2p1(n2) = 1;
            dTN3_3n1(1) = 1;
            dTN3_3p1(n3) = 1;
            dTN1_4n1(n1+1) = 1;
            dTN1_4p1(2*n1) = 1;
            dTN2_5n1(n2+1) = 1;
            dTN2_5p1(2*n2) = 1;
            dTN3_6n1(n3+1) = 1;
            dTN3_6p1(2*n3) = 1;
            dTN1_7n1(2*n1+1) = 1;
            dTN1_7p1(end-1) = 1;
            dTN2_8n1(2*n2+1) = 1;
            dTN2_8p1(end-1) = 1;
            dTN3_9n1(2*n3+1) = 1;
            dTN3_9p1(end-1) = 1;


            % building the nonlinear operator N
            % and the linear operator L

            if sqrt(kappa)*b <= 1
                N = @(v) [ D11*v - v(end).*cos(D03*v)
                    D12*v - v(end).*sin(D03*v)
                    (D01*v).*(D13*v) + v(end).*sin(D03*v) - kappa*v(end).*(D02*v).*(D01*v)
                    dT1n1*v + b
                    dT1p1*v - b
                    dT3n1*v + psib
                    dT3p1*v - psib ];

                L = @(v) [ D1, Z0, spdiags(v(end)*sin(D03*v),0,n-1,n-1)*D0, -cos(D03*v)
                    Z0, D1, spdiags(-v(end)*cos(D03*v),0,n-1,n-1)*D0, -sin(D03*v)
                    spdiags(D13*v - kappa*v(end).*(D02*v),0,n-1,n-1)*D0, spdiags(-kappa*v(end)*(D01*v),0,n-1,n-1)*D0, spdiags(D01*v,0,n-1,n-1)*D1 + spdiags(v(end)*cos(D03*v),0,n-1,n-1)*D0, sin(D03*v) - kappa*(D02*v).*(D01*v)
                    dT1n1
                    dT1p1
                    dT3n1
                    dT3p1 ];
            else

                % basic building blocks for L
                % sparse blocks created as L(row,collumn)

                L11 = spalloc((n1-1)+(n2-1)+(n3-1),n1+n2+n3,((n1-1)*n1)+((n2-1)*n2)+((n3-1)*n3));
                L11(1:n1-1,1:n1) = D1;
                L11((n1-1)+1:((n1-1)+(n2-1)),n1+1:n1+n2) = D10;
                L11((n1-1)+(n2-1)+1:((n1-1)+(n2-1)+(n3-1)),n1+n2+1:n1+n2+n3) = D100;

                L12 = spalloc((n1-1)+(n2-1)+(n3-1),n1+n2+n3,((n1-1)*n1)+((n2-1)*n2)+((n3-1)*n3));
                L12(1:n1-1,1:n1) = Z0;
                L12((n1-1)+1:((n1-1)+(n2-1)),n1+1:n1+n2) = Z00;
                L12((n1-1)+(n2-1)+1:((n1-1)+(n2-1)+(n3-1)),n1+n2+1:n1+n2+n3) = Z000;

                L13 = spalloc((n1-1)+(n2-1)+(n3-1),n1+n2+n3,((n1-1)*n1)+((n2-1)*n2)+((n3-1)*n3));
                L13(1:n1-1,1:n1) = spdiags(v(end-2)*sin(D07*v),0,n1-1,n1-1)*D0;
                L13((n1-1)+1:((n1-1)+(n2-1)),n1+1:n1+n2) = spdiags(v(end-1)*sin(D08*v),0,n2-1,n2-1)*D00;
                L13((n1-1)+(n2-1)+1:((n1-1)+(n2-1)+(n3-1)),n1+n2+1:n1+n2+n3) = spdiags(v(end)*sin(D09*v),0,n3-1,n3-1)*D000;

                L14 = spalloc((n1-1)+(n2-1)+(n3-1),3,((n1-1)+(n2-1)+(n3-1))*1);
                L14(1:n1-1,1:1) = -cos(D07*v);
                L14((n1-1)+1:((n1-1)+(n2-1)),2:2) = -cos(D08*v);
                L14((n1-1)+(n2-1)+1:((n1-1)+(n2-1)+(n3-1)),3:3) = -cos(D09*v);

                L21 = spalloc((n1-1)+(n2-1)+(n3-1),n1+n2+n3,((n1-1)*n1)+((n2-1)*n2)+((n3-1)*n3));
                L21(1:n1-1,1:n1) = Z0;
                L21((n1-1)+1:((n1-1)+(n2-1)),n1+1:n1+n2) = Z00;
                L21((n1-1)+(n2-1)+1:((n1-1)+(n2-1)+(n3-1)),n1+n2+1:n1+n2+n3) = Z000;

                L22 = spalloc((n1-1)+(n2-1)+(n3-1),n1+n2+n3,((n1-1)*n1)+((n2-1)*n2)+((n3-1)*n3));
                L22(1:n1-1,1:n1) = D1;
                L22((n1-1)+1:((n1-1)+(n2-1)),n1+1:n1+n2) = D10;
                L22((n1-1)+(n2-1)+1:((n1-1)+(n2-1)+(n3-1)),n1+n2+1:n1+n2+n3) = D100;

                L23 = spalloc((n1-1)+(n2-1)+(n3-1),n1+n2+n3,((n1-1)*n1)+((n2-1)*n2)+((n3-1)*n3));
                L23(1:n1-1,1:n1) = spdiags(-v(end-2)*cos(D07*v),0,n1-1,n1-1)*D0;
                L23((n1-1)+1:((n1-1)+(n2-1)),n1+1:n1+n2) = spdiags(-v(end-1)*cos(D08*v),0,n2-1,n2-1)*D00;
                L23((n1-1)+(n2-1)+1:((n1-1)+(n2-1)+(n3-1)),n1+n2+1:n1+n2+n3) = spdiags(-v(end)*cos(D09*v),0,n3-1,n3-1)*D000;

                L24 = spalloc((n1-1)+(n2-1)+(n3-1),3,((n1-1)+(n2-1)+(n3-1))*1);
                L24(1:n1-1,1:1) = -sin(D07*v);
                L24((n1-1)+1:((n1-1)+(n2-1)),2:2) = -sin(D08*v);
                L24((n1-1)+(n2-1)+1:((n1-1)+(n2-1)+(n3-1)),3:3) = -sin(D09*v);

                L31 = spalloc((n1-1)+(n2-1)+(n3-1),n1+n2+n3,((n1-1)*n1)+((n2-1)*n2)+((n3-1)*n3));
                L31(1:n1-1,1:n1) = spdiags(-v(end-2)*sin(D07*v)./((D01*v).^2),0,n1-1,n1-1)*D0;
                L31((n1-1)+1:((n1-1)+(n2-1)),n1+1:n1+n2) = spdiags(-v(end-1)*sin(D08*v)./((D02*v).^2),0,n2-1,n2-1)*D00;
                L31((n1-1)+(n2-1)+1:((n1-1)+(n2-1)+(n3-1)),n1+n2+1:n1+n2+n3) = spdiags(-v(end)*sin(D09*v)./((D03*v).^2),0,n3-1,n3-1)*D000;

                L32 = spalloc((n1-1)+(n2-1)+(n3-1),n1+n2+n3,((n1-1)*n1)+((n2-1)*n2)+((n3-1)*n3));
                L32(1:n1-1,1:n1) = spdiags(-kappa*v(end-2))*D0;
                L32((n1-1)+1:((n1-1)+(n2-1)),n1+1:n1+n2) = spdiags(-kappa*v(end-1))*D00;
                L32((n1-1)+(n2-1)+1:((n1-1)+(n2-1)+(n3-1)),n1+n2+1:n1+n2+n3) = spdiags(-kappa*v(end))*D000;

                L33 = spalloc((n1-1)+(n2-1)+(n3-1),n1+n2+n3,((n1-1)*n1)+((n2-1)*n2)+((n3-1)*n3));
                L33(1:n1-1,1:n1) = D1 + (spdiags(v(end-2)*cos(D07*v)./(D01*v),0,n1-1,n1-1))*D0;
                L33((n1-1)+1:((n1-1)+(n2-1)),n1+1:n1+n2) = D10 + (spdiags(v(end-1)*cos(D08*v)./(D02*v),0,n2-1,n2-1))*D00;
                L33((n1-1)+(n2-1)+1:((n1-1)+(n2-1)+(n3-1)),n1+n2+1:n1+n2+n3) = D100 + (spdiags(v(end)*cos(D09*v)./(D03*v),0,n3-1,n3-1))*D000;

                L34 = spalloc((n1-1)+(n2-1)+(n3-1),3,((n1-1)+(n2-1)+(n3-1))*1);
                L34(1:n1-1,1:1) = (sin(D07*v)./(D01*v) - kappa*D04*v);
                L34((n1-1)+1:((n1-1)+(n2-1)),2:2) = (sin(D08*v)./(D02*v) - kappa*D05*v);
                L34((n1-1)+(n2-1)+1:((n1-1)+(n2-1)+(n3-1)),3:3) = (sin(D09*v)./(D03*v) - kappa*D06*v);

                N = @(v) [ D11*v - v(end-2).*cos(D07*v)
                    D12*v - v(end-1).*cos(D08*v)
                    D13*v - v(end).*cos(D09*v)
                    D14*v - v(end-2).*sin(D07*v)
                    D15*v - v(end-1).*sin(D08*v)
                    D16*v - v(end).*sin(D09*v)
                    D17*v + v(end-2).*sin(D07*v)./(D01*v) - kappa*v(end-2).*D04*v
                    D18*v + v(end-1).*sin(D08*v)./(D02*v) - kappa*v(end-1).*D05*v
                    D19*v + v(end).*sin(D09*v)./(D03*v) - kappa*v(end).*D06*v
                    dT1n1*v + b
                    dT1p1*v + b - d
                    dT7n1*v + psib
                    dT7p1*v - dT8n1*v
                    dT2n1*v + b - d
                    dT4p1*v - dT5n1*v
                    dT2p1*v + d - b
                    dT8p1*v - dT9n1*v
                    dT3n1*v + d - b
                    dT5p1*v - dT6n1*v
                    dT3p1*v - b
                    dT9p1*v - psib ];

                L = @(v) [ L11, L12, L13, L14
                    L21, L22, L23, L24
                    L31, L32, L33, L34
                    dT1n1
                    dT1p1
                    dT7n1
                    dT7p1 - dT8n1
                    dT2n1
                    dT4p1 - dT5n1
                    dT2p1
                    dT8p1 - dT9n1
                    dT3n1
                    dT5p1 - dT6n1
                    dT3p1
                    dT9p1 ];

            end

            % initialize a counter and
            % the residual for the Newton's method loop
            kk = 1;
            res_newton = 1;

            %% Newton's method loop
            while res_newton > tol_newton

                lastwarn('', '');
                % the main steps
                %full(L(v))
                %pause
                dv = -L(v)\N(v);

                % warning syntax to catch badly scaled matrices, which
                % happened while developing the code and should no longer
                % happen
                [warnMsg, warnID] = lastwarn();
                if(isempty(warnID))
                else
                    warning('Radii and inclination angles lead to a muilt-scale problem.')
                    full(L(v));
                    pause
                    return
                    %                     % plot the current configuration if there is a problem.
                    %                     R10 = v(1:n)
                    %                     U10 = v(n+1:2*n)
                    %                     Psi10 = v(2*n+1:end-1)
                    %                     ell10 = v(end)
                    %                     figure(12)
                    %                     plot(chebfun(R10),chebfun(U10),'.-k')
                    %                     axis equal
                    %                     format longe
                    %
                    %                     pause
                end

                % Newton's step
                v = v + dv;
                %{
                plot the current configuration if there is a problem.
                tR1 = v(1:n1);
                tR2 = v(n1+1:n1+n2);
                tR3 = v(n1+n2+1:n1+n2+n3);
                tU1 = v(n1+n2+n3+1:2*n1+n2+n3);
                tU2 = v(2*n1+n2+n3+1:2*n1+2*n2+n3);
                tU3 = v(2*n1+2*n2+n3+1:2*n1+2*n2+2*n3);
                tPsi1 = v(2*n1+2*n2+2*n3+1:3*n1+2*n2+2*n3);
                tPsi2 = v(3*n1+2*n2+2*n3+1:3*n1+3*n2+2*n3);
                tPsi3 = v(3*n1+3*n2+2*n3+1:3*n1+3*n2+3*n3);
                tell01 = v(end-2);
                tell02 = v(end-1);
                tell03 = v(end);
                figure(12)
                plot(chebfun(tR1),chebfun(tU1),'.-k')
                hold on
                plot(chebfun(tR2),chebfun(tU2),'.-k')
                plot(chebfun(tR3),chebfun(tU3),'.-k')
                axis equal
                hold off

                pause
                %}

                % Build a barrier if the inclination angle leaves the
                % region that has a solution

                temp_psi = v(2*n1+2*n2+2*n3+1:3*n1+3*n2+3*n3);
                temp_psi(temp_psi > 3.5) = pi;
                temp_psi(temp_psi < -3.5) = -pi;
                v(2*n1+2*n2+2*n3+1:3*n1+3*n2+3*n3) = temp_psi;


                % computing the residual of the Newton's method
                res_newton = norm(dv,'fro')/norm(v,'fro');

                % updating the counter and if it grows above the prescribed
                % mazimum, stop this loop and report to the outer loop.
                kk = kk+1;

                if kk > max_iter_newton
                    disp('Maximum number of Newton iterations reached')
                    res_newton = res_newton
                    break
                end

            end

            %% residual and other loop conditions

            [v1, v2, v3] = extractor(v, n1, n2, n3);

            if sqrt(kappa)*b > 1

                N1 = @(v1) [ DN1_11*v1 - v1(end).*cos(DN1_07*v1)
                    DN1_14*v1 - v1(end).*sin(DN1_07*v1)
                    DN1_17*v1 + v1(end).*sin(DN1_07*v1)./(DN1_01*v1) - kappa*v1(end).*DN1_04*v1
                    dTN1_1n1*v1 + b
                    dTN1_1p1*v1 + b - d
                    dTN1_7n1*v1 + psib
                    dTN1_7p1*v1 - dTN2_8n1*v2
                    dTN1_4p1*v1 - dTN2_5n1*v2];

                N2 = @(v2) [ DN2_12*v2 - v2(end).*cos(DN2_08*v2)
                    DN2_15*v2 - v2(end).*sin(DN2_08*v2)
                    DN2_18*v2 + v2(end).*sin(DN2_08*v2)./(DN2_02*v2) - kappa*v2(end).*DN2_05*v2
                    dTN1_7p1*v1 - dTN2_8n1*v2
                    dTN2_2n1*v2 + b - d
                    dTN1_4p1*v1 - dTN2_5n1*v2
                    dTN2_2p1*v2 + d - b
                    dTN2_8p1*v2 - dTN3_9n1*v3
                    dTN2_5p1*v2 - dTN3_6n1*v3];

                N3 = @(v3) [ DN3_13*v3 - v3(end).*cos(DN3_09*v3)
                    DN3_16*v3 - v3(end).*sin(DN3_09*v3)
                    DN3_19*v3 + v3(end).*sin(DN3_09*v3)./(DN3_03*v3) - kappa*v3(end).*DN3_06*v3
                    dTN2_8p1*v2 - dTN3_9n1*v3
                    dTN3_3n1*v3 + d - b
                    dTN2_5p1*v2 - dTN3_6n1*v3
                    dTN3_3p1*v3 - b
                    dTN3_9p1*v3 - psib ];

            end

            % the relative norm of the residual

            res_bvp = norm(N(v),'fro')/norm(v,'fro');

            res_bvp1 = norm(N1(v1),'fro')/norm(v1,'fro');
            res_bvp2 = norm(N2(v2),'fro')/norm(v2,'fro');
            res_bvp3 = norm(N3(v3),'fro')/norm(v3,'fro');

            % adaptively add more Chebyshev points and resample the state
            % of the problem if the tolerance is not met.
            % Additionally, if there is too much numerical oscillation, add
            % more Chebyshev points and resample the state
            % of the problem and reset the residual.
            L21 = diffmat([2*n1 n1],0,X);
            L22 = diffmat([2*n2 n2],0,X);
            L23 = diffmat([2*n3 n3],0,X);

            c1 = length(find(diff(sign(diff((L21*v1(2*n1 + 1:3*n1)))))));
            c2 = length(find(diff(sign(diff((L22*v2(2*n2 + 1:3*n2)))))));
            c3 = length(find(diff(sign(diff((L23*v3(2*n3 + 1:3*n3)))))));

            if res_bvp1 > tol_bvp
                nn1 = n1;
                n1 = n1 + 4;
                % Sample solutions on the new grid
                S00 = diffmat([n1 nn1],0,X);
                S01 = spalloc(n1,3*nn1 + 1,(n1*nn1));
                S01(1:n1,1:nn1) = S00;
                S04 = spalloc(n1,3*nn1 + 1,(n1*nn1));
                S04(:,nn1+1:2*nn1) = S00;
                S07 = spalloc(n1,3*nn1 + 1,(n1*nn1));
                S07(:,2*nn1+1:3*nn1) = S00;
                vv1 = zeros(3*n1+1,1);
                vv1(1:n1) = S01*v1;
                vv1(n1+1:2*n1) = S04*v1;
                vv1(2*n1+1:3*n1) = S07*v1;
                vv1(end) = v1(end);
                v1 = vv1;
                [v] = compressor(v1, v2, v3);

            elseif c1 >= 2
                nn1 = n1;
                n1 = 2*(n1 - 1) -1;
                % Sample solutions on the new grid
                S001 = diffmat([n1 nn1],0,X);
                S1 = spalloc(n1,3*nn1 + 1,(n1*nn1));
                S1(:,1:nn1) = S001;
                S4 = spalloc(n1,3*nn1 + 1,(n1*nn1));
                S4(:,nn1+1:2*nn1) = S001;
                S7 = spalloc(n1,3*nn1 + 1,(n1*nn1));
                S7(:,2*nn1+1:3*nn1) = S001;
                vv1 = zeros(3*n1+1,1);
                vv1(1:n1) = S1*v1;
                vv1(n1+1:2*n1) = S4*v1;
                vv1(2*n1+1:3*n1) = S7*v1;
                vv1(end) = v1(end);
                v1 = vv1;
                res_bvp1 = 1;
                [v] = compressor(v1, v2, v3);
            end

            % if the function exceeds the maximum number of iterations,
            % break with an error statement.

            if n1 > max_iter_bvp
                disp('Maximum number of Chebyshev points reached')
                res_bvp1 = res_bvp1
                break
            end

            if res_bvp2 > tol_bvp
                nn2 = n2;
                n2 = n2 + 4;
                % Sample solutions on the new grid
                S000 = diffmat([n2 nn2],0,X);
                S02 = spalloc(n2,3*nn2 + 1,(n2*nn2));
                S02(1:n2,1:nn2) = S000;
                S05 = spalloc(n2,3*nn2 + 1,(n2*nn2));
                S05(:,nn2+1:2*nn2) = S000;
                S08 = spalloc(n2,3*nn2 + 1,(n2*nn2));
                S08(:,2*nn2+1:3*nn2) = S000;
                vv2 = zeros(3*n2+1,1);
                vv2(1:n2) = S02*v2;
                vv2(n2+1:2*n2) = S05*v2;
                vv2(2*n2+1:3*n2) = S08*v2;
                vv2(end) = v2(end);
                v2 = vv2;
                [v] = compressor(v1, v2, v3);

            elseif c2 >= 2
                nn2 = n2;
                n2 = 2*(n2 - 1) -1;
                % Sample solutions on the new grid
                S002 = diffmat([n2 nn2],0,X);
                S2 = spalloc(n2,3*nn2 + 1,(n2*nn2));
                S2(:,1:nn2) = S002;
                S5 = spalloc(n2,3*nn2 + 1,(n2*nn2));
                S5(:,nn2+1:2*nn2) = S002;
                S8 = spalloc(n2,3*nn2 + 1,(n2*nn2));
                S8(:,2*nn2+1:3*nn2) = S002;
                vv2 = zeros(3*n2+1,1);
                vv2(1:n2) = S2*v2;
                vv2(n2+1:2*n2) = S5*v2;
                vv2(2*n2+1:3*n2) = S8*v2;
                vv2(end) = v2(end);
                v2 = vv2;
                res_bvp2 = 1;
                [v] = compressor(v1, v2, v3);
            end

            % if the function exceeds the maximum number of iterations,
            % break with an error statement.

            if n2 > max_iter_bvp
                disp('Maximum number of Chebyshev points reached')
                res_bvp2 = res_bvp2
                break
            end

            if res_bvp3 > tol_bvp
                nn3 = n3;
                n3 = n3 + 4;
                % Sample solutions on the new grid
                S0000 = diffmat([n3 nn3],0,X);
                S03 = spalloc(n3,3*nn3 + 1,(n3*nn3));
                S03(1:n3,1:nn3) = S0000;
                S06 = spalloc(n3,3*nn3 + 1,(n3*nn3));
                S06(:,nn3+1:2*nn3) = S0000;
                S09 = spalloc(n3,3*nn3 + 1,(n3*nn3));
                S09(:,2*nn3+1:3*nn3) = S0000;
                vv3 = zeros(3*n3+1,1);
                vv3(1:n3) = S03*v3;
                vv3(n3+1:2*n3) = S06*v3;
                vv3(2*n3+1:3*n3) = S09*v3;
                vv3(end) = v3(end);
                v3 = vv3;
                [v] = compressor(v1, v2, v3);

            elseif  c3 >= 2
                nn3 = n3;
                n3 = 2*(n3 - 1) -1;
                % Sample solutions on the new grid
                S003 = diffmat([n3 nn3],0,X);
                S3 = spalloc(n3,3*nn3 + 1,(n3*nn3));
                S3(:,1:nn3) = S003;
                S6 = spalloc(n3,3*nn3 + 1,(n3*nn3));
                S6(:,nn3+1:2*nn3) = S003;
                S9 = spalloc(n3,3*nn3 + 1,(n3*nn3));
                S9(:,2*nn3+1:3*nn3) = S003;
                vv3 = zeros(3*n3+1,1);
                vv3(1:n3) = S3*v3;
                vv3(n3+1:2*n3) = S6*v3;
                vv3(2*n3+1:3*n3) = S9*v3;
                vv3(end) = v3(end);
                v3 = vv3;
                res_bvp3 = 1;
                [v] = compressor(v1, v2, v3);
            end

            % if the function exceeds the maximum number of iterations,
            % break with an error statement.

            if n3 > max_iter_bvp
                disp('Maximum number of Chebyshev points reached')
                res_bvp3 = res_bvp3
                break
            end
        end
    end

%% the main computational engine of the file
    function [v, n] = cheb_engine_small(v, n)

        % intialize the residual
        res_bvp = 1;

        while res_bvp > tol_bvp

            % initialize the differential operator components
            %
            % D0 and D1 are spectral operators corresponding to a
            % downsampling matrix and a differentiation matrix,
            % respectively
            D0 = diffmat([n-1 n],0,X);
            D1 = diffmat([n-1 n],1,X);
            Z0 = sparse(n-1,n);
            D01 = spalloc(n-1, 3*n + 1, n*(n-1));
            D02 = D01;
            D03 = D01;
            D11 = D01;
            D12 = D01;
            D13 = D01;
            D01(1:n-1, 1:n) = D0;
            D11(1:n-1, 1:n) = D1;
            D02(1:n-1, n+1:2*n) = D0;
            D12(1:n-1, n+1:2*n) = D1;
            D03(1:n-1, 2*n+1:3*n) = D0;
            D13(1:n-1, 2*n+1:3*n) = D1;

            % Evaluate the computational vector to check the boundary
            % conditions
            dT1n1 = sparse(1,3*n+1);
            dT1p1 = dT1n1;
            dT3n1 = dT1n1;
            dT3p1 = dT1n1;
            dT1n1(1) = 1;
            dT1p1(n) = 1;
            dT3n1(2*n+1) = 1;
            dT3p1(end -1) = 1;

            % building the nonlinear operator N
            % and the linear operator L

            if sqrt(kappa)*b <= 1
                N = @(v) [ D11*v - v(end).*cos(D03*v)
                    D12*v - v(end).*sin(D03*v)
                    (D01*v).*(D13*v) + v(end).*sin(D03*v) - kappa*v(end).*(D02*v).*(D01*v)
                    dT1n1*v + b
                    dT1p1*v - b
                    dT3n1*v + psib
                    dT3p1*v - psib ];

                L = @(v) [ D1, Z0, spdiags(v(end)*sin(D03*v),0,n-1,n-1)*D0, -cos(D03*v)
                    Z0, D1, spdiags(-v(end)*cos(D03*v),0,n-1,n-1)*D0, -sin(D03*v)
                    spdiags(D13*v - kappa*v(end).*(D02*v),0,n-1,n-1)*D0, spdiags(-kappa*v(end)*(D01*v),0,n-1,n-1)*D0, spdiags(D01*v,0,n-1,n-1)*D1 + spdiags(v(end)*cos(D03*v),0,n-1,n-1)*D0, sin(D03*v) - kappa*(D02*v).*(D01*v)
                    dT1n1
                    dT1p1
                    dT3n1
                    dT3p1 ];
            else

                N = @(v) [ D11*v - v(end).*cos(D03*v)
                    D12*v - v(end).*sin(D03*v)
                    D13*v + v(end).*sin(D03*v)./(D01*v) - kappa*v(end).*D02*v
                    dT1n1*v + b
                    dT1p1*v - b
                    dT3n1*v + psib
                    dT3p1*v - psib ];

                L = @(v) [ D1, Z0, spdiags(v(end)*sin(D03*v),0,n-1,n-1)*D0, -cos(D03*v)
                    Z0, D1, spdiags(-v(end)*cos(D03*v),0,n-1,n-1)*D0, -sin(D03*v)
                    spdiags(-v(end)*sin(D03*v)./((D01*v).^2),0,n-1,n-1)*D0, -kappa*v(end)*D0, D1 + (spdiags(v(end)*cos(D03*v)./(D01*v),0,n-1,n-1))*D0, sin(D03*v)./(D01*v) - kappa*D02*v
                    dT1n1
                    dT1p1
                    dT3n1
                    dT3p1 ];
            end

            % initialize a counter and
            % the residual for the Newton's method loop
            kk = 1;
            res_newton = 1;

            %% Newton's method loop
            while res_newton > tol_newton

                lastwarn('', '');
                % the main steps
                dv = -L(v)\N(v);

                % warning syntax to catch badly scaled matrices, which
                % happened while developing the code and should no longer
                % happen
                [warnMsg, warnID] = lastwarn();
                if(isempty(warnID))
                else
                    warning('Radii and inclination angles lead to a muilt-scale problem.')
                    return
                    %                     % plot the current configuration if there is a problem.
                    %                     R10 = v(1:n)
                    %                     U10 = v(n+1:2*n)
                    %                     Psi10 = v(2*n+1:end-1)
                    %                     ell10 = v(end)
                    %                     figure(12)
                    %                     plot(chebfun(R10),chebfun(U10),'.-k')
                    %                     axis equal
                    %                     format longe
                    %
                    %                     pause
                end

                % Newton's step
                v = v + dv;

                % Build a barrier if the inclination angle leaves the
                % region that has a solution

                temp_psi = v(2*n+1:3*n);
                temp_psi(temp_psi > 3.5) = pi;
                temp_psi(temp_psi < -3.5) = -pi;
                v(2*n+1:3*n) = temp_psi;


                % computing the residual of the Newton's method
                res_newton = norm(dv,'fro')/norm(v,'fro');

                % updating the counter and if it grows above the prescribed
                % mazimum, stop this loop and report to the outer loop.
                kk = kk+1;

                if kk > max_iter_newton
                    disp('Maximum number of Newton iterations reached')
                    res_newton = res_newton
                    break
                end

            end

            %% residual and other loop conditions

            % the relative norm of the residual
            res_bvp = norm(N(v),'fro')/norm(v,'fro');

            % adaptively add more Chebyshev points and resample the state
            % of the problem if the tolerance is not met.
            % Additionally, if there is too much numerical oscillation, add
            % more Chebyshev points and resample the state
            % of the problem and reset the residual.
            S2 = diffmat([2*n n],0,X);
            if res_bvp > tol_bvp
                nn = n;
                n = n + 4;
                % Sample solutions on the new grid
                S0 = diffmat([n nn],0,X);
                S01 = spalloc(n, 3*nn + 1,n*nn);
                S01(:,1:nn) = S0;
                S02 = spalloc(n, 3*nn + 1,n*nn);
                S02(:,nn+1:2*nn) = S0;
                S03 = spalloc(n, 3*nn + 1,n*nn);
                S03(:,2*nn+1:3*nn) = S0;
                vv = zeros(3*n+1,1);
                vv(1:n) = S01*v;
                vv(n+1:2*n) = S02*v;
                vv(2*n+1:3*n) = S03*v;
                vv(end) = v(end);
                v = vv;

            elseif length(find(diff(sign(diff((S2*v(2*n+1:3*n))))))) >= 2
                nn = n;
                n = 2*(n - 1) -1;
                % Sample solutions on the new grid
                S0 = diffmat([n nn],0,X);
                S01 = spalloc(n, 3*nn + 1,n*nn);
                S01(:,1:nn) = S0;
                S02 = spalloc(n, 3*nn + 1,n*nn);
                S02(:,nn+1:2*nn) = S0;
                S03 = spalloc(n, 3*nn + 1,n*nn);
                S03(:,2*nn+1:3*nn) = S0;
                vv = zeros(3*n+1,1);
                vv(1:n) = S01*v;
                vv(n+1:2*n) = S02*v;
                vv(2*n+1:3*n) = S03*v;
                vv(end) = v(end);
                v = vv;

                res_bvp = 1;
            else
                break
            end

            % if the function exceeds the maximum number of iterations,
            % break with an error statement.
            if n > max_iter_bvp
                disp('Maximum number of Chebyshev points reached')
                res_bvp = res_bvp
                break
            end
        end
    end

%% The function for exctracting v1, v2, v3

    function [v1, v2, v3] = extractor(v, n1, n2, n3)

        v1 = zeros(3*n1+1,1);
        v1(1:n1) = v(1:n1);
        v1(n1+1:2*n1) = v(n1+n2+n3+1:2*n1+n2+n3);
        v1(2*n1+1:3*n1) = v(2*n1+2*n2+2*n3+1:3*n1+2*n2+2*n3);
        v1(end) = v(end-2);

        v2 = zeros(3*n2+1,1);
        v2(1:n2) = v(n1+1:n1+n2);
        v2(n2+1:2*n2) = v(2*n1+n2+n3+1:2*n1+2*n2+n3);
        v2(2*n2+1:3*n2) = v(3*n1+2*n2+2*n3+1:3*n1+3*n2+2*n3);
        v2(end) = v(end-1);

        v3 = zeros(3*n3+1,1);
        v3(1:n3) = v(n1+n2+1:n1+n2+n3);
        v3(n3+1:2*n3) = v(2*n1+2*n2+n3+1:2*n1+2*n2+2*n3);
        v3(2*n3+1:3*n3) = v(3*n1+3*n2+2*n3+1:3*n1+3*n2+3*n3);
        v3(end) = v(end);

    end

%% The function for compressing v1, v2, v3

    function [v] = compressor(v1, v2, v3)

        cn1 = (length(v1) - 1)/3;
        cn2 = (length(v2) - 1)/3;
        cn3 = (length(v3) - 1)/3;
        v = zeros(3*(cn1+cn2+cn3+1),1);
        v(1:cn1) = v1(1:cn1);
        v(cn1+1:cn1+cn2) = v2(1:cn2);
        v(cn1+cn2+1:cn1+cn2+cn3) = v3(1:cn3);
        v(cn1+cn2+cn3+1:2*cn1+cn2+cn3) = v1(cn1+1:2*cn1);
        v(2*cn1+cn2+cn3+1:2*cn1+2*cn2+cn3) = v2(cn2+1:2*cn2);
        v(2*cn1+2*cn2+cn3+1:2*cn1+2*cn2+2*cn3) = v3(cn3+1:2*cn3);
        v(2*cn1+2*cn2+2*cn3+1:3*cn1+2*cn2+2*cn3) = v1(2*cn1+1:3*cn1);
        v(3*cn1+2*cn2+2*cn3+1:3*cn1+3*cn2+2*cn3) = v2(2*cn2+1:3*cn2);
        v(3*cn1+3*cn2+2*cn3+1:3*cn1+3*cn2+3*cn3) = v3(2*cn3+1:3*cn3);
        v(end-2) = v1(end);
        v(end-1) = v2(end);
        v(end) = v3(end);




    end


%% Assigning the output into variables that have a physical meaning
if sqrt(kappa)*b >= 1
    R1 = v(1:n1);
    R2 = v(n1+1:n1+n2);
    R3 = v(n1+n2+1:n1+n2+n3);
    U1 = v(n1+n2+n3+1:2*n1+n2+n3);
    U2 = v(2*n1+n2+n3+1:2*n1+2*n2+n3);
    U3 = v(2*n1+2*n2+n3+1:2*n1+2*n2+2*n3);
    Psi1 = v(2*n1+2*n2+2*n3+1:3*n1+2*n2+2*n3);
    Psi2 = v(3*n1+2*n2+2*n3+1:3*n1+3*n2+2*n3);
    Psi3 = v(3*n1+3*n2+2*n3+1:3*n1+3*n2+3*n3);
    ell01 = v(end-2);
    ell02 = v(end-1);
    ell03 = v(end);

    %% Plotting
    % delete or comment this block if the function is to be used inside some
    % other algortihm.
    figure(2)
    plot(chebfun(R1),chebfun(U1),'.-k')
    hold on
    plot(chebfun(R2),chebfun(U2),'.-k')
    plot(chebfun(R3),chebfun(U3),'.-k')
    axis equal
    plot([-b;-b],ylim,'-.k')
    plot([b;b],ylim,'-.k')
    hold off

else
    R = v(1:n);
    U = v(n+1:2*n);
    Psi = v(2*n+1:end-1);
    ell = v(end);

    %% Plotting
    % delete or comment this block if the function is to be used inside some
    % other algortihm.
    figure(1)
    plot(chebfun(R),chebfun(U),'.-k')
    hold on
    axis equal
    plot([-b;-b],ylim,'-.k')
    plot([b;b],ylim,'-.k')
end
end