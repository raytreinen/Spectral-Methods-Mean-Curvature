function [] = symmetric_capillary_annular_2z_to_3z()
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
% function [R, U, Psi ,ell, n, res_bvp] = symmetric_capillary_annular(a,b,psia_actual,psib_actual,kappa)
% and delete (or comment) the physical parameters block as well as the figure plotting
% at the end of the program.

%% physical parameters
% this block can be commented out if the program is used inside of a larger
% program and these values are passed as inputs to this program.
kappa = 1;
a = .05;
b = 5;
d = .2;
psia_actual = -31*pi/32;
psib_actual = -16*pi/16;
psic = -pi/2;

%tic
%% Computational parameters
% computational domain
X = [-1,1];

% maximum loop counts and the tolerances for each loop
max_iter_newton = 100000;
max_iter_bvp = 10000;
tol_newton = 1e-13;
tol_bvp = 1e-12;

% initialize the number of points.
k = 7;
n = 2*k + 1;
n1 = 2*k + 1;
n2 = 2*k + 1;
n3 = k + 4;

% number of breaks for angle with abs value greater than pi/2
m = 21;

%% Initial guesses

% assign to psia and psib the actual values if they are less than pi/2 in
% magnitude, and if not, then replace them with +- pi/2
psia = sign(psia_actual)*min(abs(psia_actual), pi/2);
psib = sign(psib_actual)*min(abs(psib_actual), pi/2);


% compute the Chebyshev points (called tau in the paper)
s1 = chebpts(n1);
s3 = chebpts(n3);

% compute initial guesses for psi, R and U.  These may be replaced below.
Psi01 = @(s1) (1 - s1)*psia/2;
Psi03 = @(s3) (1 + s3)*psib/2;

if (psia*psib <= 0) && (psib >= 0)
    [R, U, Psi, ell] = symmetric_capillary_annular_interpolation_func(a,b,d,-pi/2,pi/2,kappa,k);
elseif (psia*psib <= 0) && (psib < 0)
    [R, U, Psi, ell] = symmetric_capillary_annular_interpolation_func(a,b,d,pi/2,-pi/2,kappa,k);
elseif (psia*psib >= 0) && (psib >= 0)
    [R, U, Psi, ell] = symmetric_capillary_annular_interpolation_func(a,b,d,pi/2,pi/2,kappa,k)
else
    [R, U, Psi, ell] = symmetric_capillary_annular_interpolation_func(a,b,d,-pi/2,-pi/2,kappa,k);
end

s = chebpts(length(R));
R0 = chebfun(R);
U0 = chebfun(U);
Psi0 = chebfun(Psi);

%% Interpolation

% Computing Zone 1
ss = interp1(R,s,a+d);
si_1 = chebpts(n1,[-1,1]);
Ri_1 = R0((((ss+1)/2)*(si_1+1))-1);
Ui_1 = U0((((ss+1)/2)*(si_1+1))-1);
Pi_1 = Psi0((((ss+1)/2)*(si_1+1))-1);

ell0_int_1 = (ell/2)*(ss+1);

% Computing Zone 2
si_2 = chebpts(n3,[-1,1]);
Ri_2 = R0((((1-ss)/2)*(si_2-1))+1);
Ui_2 = U0((((1-ss)/2)*(si_2-1))+1);
Pi_2 = Psi0((((1-ss)/2)*(si_2-1))+1);

ell0_int_2 = ell - ell0_int_1;

% Interpolated 'v'
v = [Ri_1; Ri_2; Ui_1; Ui_2; Pi_1; Pi_2; ell0_int_1; ell0_int_2];

%% solving the problem if it is a graph over a base domain and configuring the problem if it is not

% tracking the newton iterations
newt_numb = 0;

% solve the problem if the solution is a graph over a base domain
[v, n1, n3] = cheb_engine_med(v, n1, n3);

% solve the problem if the solution is not a graph over a base domain
if ((abs(psib_actual) > pi/2)||(abs(psia_actual) > pi/2))
    psia_vec = linspace(sign(psia_actual)*pi/2,psia_actual,m)';
    psia_vec = psia_vec(2:end);
    psib_vec = linspace(sign(psib_actual)*pi/2,psib_actual,m)';
    psib_vec = psib_vec(2:end);

    for i = 1:m-1
        if i < m-2
            psia = psia_vec(i);
            psib = psib_vec(i);
            [v, n1, n3] = cheb_engine_med(v, n1, n3);
        elseif i == m-2
            psia = psia_vec(i);
            psib = psib_vec(i);
            [v, n1, n3] = cheb_engine_med(v, n1, n3);
        else
            % Break zone 1 into two zones at psic
            R1 = v(1:n1);
            R3 = v(n1+1:n1+n3);
            U1 = v(n1+n3+1:2*n1+n3);
            U3 = v(2*n1+n3+1:2*n1+2*n3);
            Psi1 = v(2*n1+2*n3+1:3*n1+2*n3);
            Psi3 = v(3*n1+2*n3+1:3*n1+3*n3);
            ell01 = v(end-1);
            ell03 = v(end);
            s1 = chebpts(n1);
            R1M = chebfun(R1);
            U1M = chebfun(U1);
            Psi1M = chebfun(Psi1);

            % points for zone within zone
            n2 = n1;
            n1 = k;

            % Interpolation
            sss = interp1(Psi1,s1,psic);

            % Computing Zone 1a
            si_1a = chebpts(n1,[-1,1]);
            Ri_1a = R1M((((sss+1)/2)*(si_1a+1))-1);
            Ui_1a = U1M((((sss+1)/2)*(si_1a+1))-1);
            Pi_1a = Psi1M((((sss+1)/2)*(si_1a+1))-1);

            ell0_int_1a = (ell01/2)*(sss+1);

            % Computing Zone 1b
            si_1b = chebpts(n2,[-1,1]);
            Ri_1b = R1M((((1-sss)/2)*(si_1b-1))+1);
            Ui_1b = U1M((((1-sss)/2)*(si_1b-1))+1);
            Pi_1b = Psi1M((((1-sss)/2)*(si_1b-1))+1);

            ell0_int_1b = ell01 - ell0_int_1a;

            % Interpolated 'v' (zone 11, zone 12, zone 3)
            i;
            v = [Ri_1a; Ri_1b; R3; Ui_1a; Ui_1b; U3; Pi_1a; Pi_1b; Psi3; ell0_int_1a; ell0_int_1b; ell03];
            % pause

            % Continuing to rachet
            psia = psia_vec(i);
            psib = psib_vec(i);
            [v, n1, n2, n3] = cheb_engine_large(v, n1, n2, n3);
        end

    end
end


%% the main computational engine of the file
    function [v, n1, n3] = cheb_engine_med(v, n1, n3)

        % intialize the residual
        res_bvp1 = 1;
        res_bvp3 = 1;

        while (res_bvp1 > tol_bvp)||(res_bvp3 > tol_bvp)

            % initialize the differential operator components
            %
            % D0 and D1 are spectral operators corresponding to a
            % downsampling matrix and a differentiation matrix,
            % respectively
            D0 = diffmat([n1-1 n1],0,X);
            D1 = diffmat([n1-1 n1],1,X);
            Z0 = sparse(n1-1,n1);
            D000 = diffmat([n3-1 n3],0,X);
            D100 = diffmat([n3-1 n3],1,X);
            Z000 = sparse(n3-1,n3);
            D01 = spalloc(n1-1, 3*(n1 + n3)+2, n1*(n1-1));
            D03 = spalloc(n3-1, 3*(n1 + n3)+2, n3*(n3-1));
            D04 = D01;
            D06 = D03;
            D07 = D01;
            D09 = D03;
            D11 = D01;
            D13 = D03;
            D14 = D01;
            D16 = D03;
            D17 = D01;
            D19 = D03;
            D01(1:n1-1, 1:n1) = D0;
            D11(1:n1-1, 1:n1) = D1;
            D03(1:n3-1, n1+1:n1+n3) = D000;
            D13(1:n3-1, n1+1:n1+n3) = D100;
            D04(1:n1-1, n1+n3+1:2*n1+n3) = D0;
            D14(1:n1-1, n1+n3+1:2*n1+n3) = D1;
            D06(1:n3-1, 2*n1+n3+1:2*n1+2*n3) = D000;
            D16(1:n3-1, 2*n1+n3+1:2*n1+2*n3) = D100;
            D07(1:n1-1, 2*n1+2*n3+1:3*n1+2*n3) = D0;
            D17(1:n1-1, 2*n1+2*n3+1:3*n1+2*n3) = D1;
            D09(1:n3-1, 3*n1+2*n3+1:3*n1+3*n3) = D000;
            D19(1:n3-1, 3*n1+2*n3+1:3*n1+3*n3) = D100;

            DN1_01 = spalloc(n1-1, 3*n1+1, n1*(n1-1));
            DN3_03 = spalloc(n3-1, 3*n3+1, n3*(n3-1));
            DN1_04 = DN1_01;
            DN3_06 = DN3_03;
            DN1_07 = DN1_01;
            DN3_09 = DN3_03;
            DN1_11 = DN1_01;
            DN3_13 = DN3_03;
            DN1_14 = DN1_01;
            DN3_16 = DN3_03;
            DN1_17 = DN1_01;
            DN3_19 = DN3_03;
            DN1_01(1:n1-1, 1:n1) = D0;
            DN1_11(1:n1-1, 1:n1) = D1;
            DN3_03(1:n3-1, 1:n3) = D000;
            DN3_13(1:n3-1, 1:n3) = D100;
            DN1_04(1:n1-1, n1+1:2*n1) = D0;
            DN1_14(1:n1-1, n1+1:2*n1) = D1;
            DN3_06(1:n3-1, n3+1:2*n3) = D000;
            DN3_16(1:n3-1, n3+1:2*n3) = D100;
            DN1_07(1:n1-1, 2*n1+1:3*n1) = D0;
            DN1_17(1:n1-1, 2*n1+1:3*n1) = D1;
            DN3_09(1:n3-1, 2*n3+1:3*n3) = D000;
            DN3_19(1:n3-1, 2*n3+1:3*n3) = D100;


            % Evaluate the computational vector to check the boundary
            % conditions
            dT1n1 = sparse(1,3*(n1 + n3)+2);
            dT1p1 = dT1n1;
            dT3n1 = dT1n1;
            dT3p1 = dT1n1;
            dT4n1 = dT1n1;
            dT4p1 = dT1n1;
            dT6n1 = dT1n1;
            dT6p1 = dT1n1;
            dT7n1 = dT1n1;
            dT7p1 = dT1n1;
            dT9n1 = dT1n1;
            dT9p1 = dT1n1;
            dT1n1(1) = 1;
            dT1p1(n1) = 1;
            dT3n1(n1+1) = 1;
            dT3p1(n1+n3) = 1;
            dT4n1(n1+n3+1) = 1;
            dT4p1(2*n1+n3) = 1;
            dT6n1(2*n1+n3+1) = 1;
            dT6p1(2*n1+2*n3) = 1;
            dT7n1(2*n1+2*n3+1) = 1;
            dT7p1(3*n1+2*n3) = 1;
            dT9n1(3*n1+2*n3+1) = 1;
            dT9p1(end-2) = 1;

            dTN1_1n1 = sparse(1,3*n1+1);
            dTN1_1p1 = dTN1_1n1;
            dTN3_3n1 = sparse(1,3*n3+1);
            dTN3_3p1 = dTN3_3n1;
            dTN1_4n1 = dTN1_1n1;
            dTN1_4p1 = dTN1_1n1;
            dTN3_6n1 = dTN3_3n1;
            dTN3_6p1 = dTN3_3n1;
            dTN1_7n1 = dTN1_1n1;
            dTN1_7p1 = dTN1_1n1;
            dTN3_9n1 = dTN3_3n1;
            dTN3_9p1 = dTN3_3n1;
            dTN1_1n1(1) = 1;
            dTN1_1p1(n1) = 1;
            dTN3_3n1(1) = 1;
            dTN3_3p1(n3) = 1;
            dTN1_4n1(n1+1) = 1;
            dTN1_4p1(2*n1) = 1;
            dTN3_6n1(n3+1) = 1;
            dTN3_6p1(2*n3) = 1;
            dTN1_7n1(2*n1+1) = 1;
            dTN1_7p1(end-1) = 1;
            dTN3_9n1(2*n3+1) = 1;
            dTN3_9p1(end-1) = 1;


            % building the nonlinear operator N
            % and the linear operator L

            if sqrt(kappa)*b <= 3

                % basic building blocks for L
                % sparse blocks created as L(row,collumn)

                L11 = spalloc((n1-1)+(n3-1),n1+n3,((n1-1)*n1)+((n3-1)*n3));
                L11(1:n1-1,1:n1) = D1;
                L11((n1-1)+1:((n1-1)+(n3-1)),n1+1:n1+n3) = D100;

                L12 = spalloc((n1-1)+(n3-1),n1+n3,((n1-1)*n1)+((n3-1)*n3));
                L12(1:n1-1,1:n1) = Z0;
                L12((n1-1)+1:((n1-1)+(n3-1)),n1+1:n1+n3) = Z000;

                L13 = spalloc((n1-1)+(n3-1),n1+n3,((n1-1)*n1)+((n3-1)*n3));
                L13(1:n1-1,1:n1) = spdiags(v(end-1)*sin(D07*v),0,n1-1,n1-1)*D0;
                L13((n1-1)+1:((n1-1)+(n3-1)),n1+1:n1+n3) = spdiags(v(end)*sin(D09*v),0,n3-1,n3-1)*D000;

                L14 = spalloc((n1-1)+(n3-1),2,((n1-1)+(n3-1))*1);
                L14(1:n1-1,1:1) = -cos(D07*v);
                L14((n1-1)+1:((n1-1)+(n3-1)),2:2) = -cos(D09*v);

                L21 = spalloc((n1-1)+(n3-1),n1+n3,((n1-1)*n1)+((n3-1)*n3));
                L21(1:n1-1,1:n1) = Z0;
                L21((n1-1)+1:((n1-1)+(n3-1)),n1+1:n1+n3) = Z000;

                L22 = spalloc((n1-1)+(n3-1),n1+n3,((n1-1)*n1)+((n3-1)*n3));
                L22(1:n1-1,1:n1) = D1;
                L22((n1-1)+1:((n1-1)+(n3-1)),n1+1:n1+n3) = D100;

                L23 = spalloc((n1-1)+(n3-1),n1+n3,((n1-1)*n1)+((n3-1)*n3));
                L23(1:n1-1,1:n1) = spdiags(-v(end-1)*cos(D07*v),0,n1-1,n1-1)*D0;
                L23((n1-1)+1:((n1-1)+(n3-1)),n1+1:n1+n3) = spdiags(-v(end)*cos(D09*v),0,n3-1,n3-1)*D000;

                L24 = spalloc((n1-1)+(n3-1),2,((n1-1)+(n3-1))*1);
                L24(1:n1-1,1:1) = -sin(D07*v);
                L24((n1-1)+1:((n1-1)+(n3-1)),2:2) = -sin(D09*v);

                L31 = spalloc((n1-1)+(n3-1),n1+n3,((n1-1)*n1)+((n3-1)*n3));
                L31(1:n1-1,1:n1) = spdiags(D17*v - kappa*v(end-1).*(D04*v),0,n1-1,n1-1)*D0;
                L31((n1-1)+1:((n1-1)+(n3-1)),n1+1:n1+n3) = spdiags(D19*v - kappa*v(end).*(D06*v),0,n3-1,n3-1)*D000;

                L32 = spalloc((n1-1)+(n3-1),n1+n3,((n1-1)*n1)+((n3-1)*n3));
                L32(1:n1-1,1:n1) = spdiags(-kappa*v(end-1)*(D01*v),0,n1-1,n1-1)*D0;
                L32((n1-1)+1:((n1-1)+(n3-1)),n1+1:n1+n3) = spdiags(-kappa*v(end)*(D03*v),0,n3-1,n3-1)*D000;

                L33 = spalloc((n1-1)+(n3-1),n1+n3,((n1-1)*n1)+((n3-1)*n3));
                L33(1:n1-1,1:n1) = spdiags(D01*v,0,n1-1,n1-1)*D1 + spdiags(v(end-1)*cos(D07*v),0,n1-1,n1-1)*D0;
                L33((n1-1)+1:((n1-1)+(n3-1)),n1+1:n1+n3) = spdiags(D03*v,0,n3-1,n3-1)*D100 + spdiags(v(end)*cos(D09*v),0,n3-1,n3-1)*D000;

                L34 = spalloc((n1-1)+(n3-1),2,((n1-1)+(n3-1))*1);
                L34(1:n1-1,1:1) = sin(D07*v) - kappa*(D04*v).*(D01*v);
                L34((n1-1)+1:((n1-1)+(n3-1)),2:2) = sin(D09*v) - kappa*(D06*v).*(D03*v);

                N = @(v) [ D11*v - v(end-1).*cos(D07*v)
                    D13*v - v(end).*cos(D09*v)
                    D14*v - v(end-1).*sin(D07*v)
                    D16*v - v(end).*sin(D09*v)
                    (D01*v).*(D17*v) + v(end-1).*sin(D07*v) - kappa*v(end-1).*(D04*v).*(D01*v)
                    (D03*v).*(D19*v) + v(end).*sin(D09*v) - kappa*v(end).*(D06*v).*(D03*v)
                    dT1n1*v - a
                    dT1p1*v - (a + d)
                    dT3n1*v - (a + d)
                    dT3p1*v - b
                    dT4p1*v - dT6n1*v
                    dT7n1*v - psia
                    dT7p1*v - dT9n1*v
                    dT9p1*v - psib ];

                L = @(v) [ L11, L12, L13, L14
                    L21, L22, L23, L24
                    L31, L32, L33, L34
                    dT1n1
                    dT1p1
                    dT3n1
                    dT3p1
                    dT4p1 - dT6n1
                    dT7n1
                    dT7p1 - dT9n1
                    dT9p1];

            else

                % basic building blocks for L
                % sparse blocks created as L(row,collumn)

                L11 = spalloc((n1-1)+(n3-1),n1+n3,((n1-1)*n1)+((n3-1)*n3));
                L11(1:n1-1,1:n1) = D1;
                L11((n1-1)+1:((n1-1)+(n3-1)),n1+1:n1+n3) = D100;

                L12 = spalloc((n1-1)+(n3-1),n1+n3,((n1-1)*n1)+((n3-1)*n3));
                L12(1:n1-1,1:n1) = Z0;
                L12((n1-1)+1:((n1-1)+(n3-1)),n1+1:n1+n3) = Z000;

                L13 = spalloc((n1-1)+(n3-1),n1+n3,((n1-1)*n1)+((n3-1)*n3));
                L13(1:n1-1,1:n1) = spdiags(v(end-1)*sin(D07*v),0,n1-1,n1-1)*D0;
                L13((n1-1)+1:((n1-1)+(n3-1)),n1+1:n1+n3) = spdiags(v(end)*sin(D09*v),0,n3-1,n3-1)*D000;

                L14 = spalloc((n1-1)+(n3-1),2,((n1-1)+(n3-1))*1);
                L14(1:n1-1,1:1) = -cos(D07*v);
                L14((n1-1)+1:((n1-1)+(n3-1)),2:2) = -cos(D09*v);

                L21 = spalloc((n1-1)+(n3-1),n1+n3,((n1-1)*n1)+((n3-1)*n3));
                L21(1:n1-1,1:n1) = Z0;
                L21((n1-1)+1:((n1-1)+(n3-1)),n1+1:n1+n3) = Z000;

                L22 = spalloc((n1-1)+(n3-1),n1+n3,((n1-1)*n1)+((n3-1)*n3));
                L22(1:n1-1,1:n1) = D1;
                L22((n1-1)+1:((n1-1)+(n3-1)),n1+1:n1+n3) = D100;

                L23 = spalloc((n1-1)+(n3-1),n1+n3,((n1-1)*n1)+((n3-1)*n3));
                L23(1:n1-1,1:n1) = spdiags(-v(end-1)*cos(D07*v),0,n1-1,n1-1)*D0;
                L23((n1-1)+1:((n1-1)+(n3-1)),n1+1:n1+n3) = spdiags(-v(end)*cos(D09*v),0,n3-1,n3-1)*D000;

                L24 = spalloc((n1-1)+(n3-1),2,((n1-1)+(n3-1))*1);
                L24(1:n1-1,1:1) = -sin(D07*v);
                L24((n1-1)+1:((n1-1)+(n3-1)),2:2) = -sin(D09*v);

                L31 = spalloc((n1-1)+(n3-1),n1+n3,((n1-1)*n1)+((n3-1)*n3));
                L31(1:n1-1,1:n1) = spdiags(D17*v - kappa*v(end-1).*(D04*v),0,n1-1,n1-1)*D0;
                L31((n1-1)+1:((n1-1)+(n3-1)),n1+1:n1+n3) = spdiags(-v(end)*sin(D09*v)./((D03*v).^2),0,n3-1,n3-1)*D000;

                L32 = spalloc((n1-1)+(n3-1),n1+n3,((n1-1)*n1)+((n3-1)*n3));
                L32(1:n1-1,1:n1) = spdiags(-kappa*v(end-1)*(D01*v),0,n1-1,n1-1)*D0;
                L32((n1-1)+1:((n1-1)+(n3-1)),n1+1:n1+n3) = spdiags(-kappa*v(end))*D000;

                L33 = spalloc((n1-1)+(n3-1),n1+n3,((n1-1)*n1)+((n3-1)*n3));
                L33(1:n1-1,1:n1) = spdiags(D01*v,0,n1-1,n1-1)*D1 + spdiags(v(end-1)*cos(D07*v),0,n1-1,n1-1)*D0;
                L33((n1-1)+1:((n1-1)+(n3-1)),n1+1:n1+n3) = D100 + (spdiags(v(end)*cos(D09*v)./(D03*v),0,n3-1,n3-1))*D000;

                L34 = spalloc((n1-1)+(n3-1),2,((n1-1)+(n3-1))*1);
                L34(1:n1-1,1:1) = sin(D07*v) - kappa*(D04*v).*(D01*v);
                L34((n1-1)+1:((n1-1)+(n3-1)),2:2) = (sin(D09*v)./(D03*v) - kappa*D06*v);

                N = @(v) [ D11*v - v(end-1).*cos(D07*v)
                    D13*v - v(end).*cos(D09*v)
                    D14*v - v(end-1).*sin(D07*v)
                    D16*v - v(end).*sin(D09*v)
                    (D01*v).*(D17*v) + v(end-1).*sin(D07*v) - kappa*v(end-1).*(D04*v).*(D01*v)
                    D19*v + v(end).*sin(D09*v)./(D03*v) - kappa*v(end).*D06*v
                    dT1n1*v - a
                    dT1p1*v - (a + d)
                    dT3n1*v - (a + d)
                    dT3p1*v - b
                    dT4p1*v - dT6n1*v
                    dT7n1*v - psia
                    dT7p1*v - dT9n1*v
                    dT9p1*v - psib ];

                L = @(v) [ L11, L12, L13, L14
                    L21, L22, L23, L24
                    L31, L32, L33, L34
                    dT1n1
                    dT1p1
                    dT3n1
                    dT3p1
                    dT4p1 - dT6n1
                    dT7n1
                    dT7p1 - dT9n1
                    dT9p1];

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

                % Build a barrier if the inclination angle leaves the
                % region that has a solution

                temp_psi = v(2*n1+2*n3+1:3*n1+3*n3);
                temp_psi(temp_psi > 3.5) = pi;
                temp_psi(temp_psi < -3.5) = -pi;
                v(2*n1+2*n3+1:3*n1+3*n3) = temp_psi;


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

            [v1, v3] = extractor_med(v, n1, n3);

            if sqrt(kappa)*b > 3

                N1 = @(v1) [ DN1_11*v1 - v1(end).*cos(DN1_07*v1)
                    DN1_14*v1 - v1(end).*sin(DN1_07*v1)
                    (DN1_01*v1).*(DN1_17*v1) + v1(end).*sin(DN1_07*v1) - kappa*v1(end).*(DN1_04*v1).*(DN1_01*v1)
                    dTN1_1n1*v1 - a
                    dTN1_1p1*v1 - (a+d)
                    dTN1_7n1*v1 - psia
                    dTN1_7p1*v1 - dTN3_9n1*v3
                    dTN1_4p1*v1 - dTN3_6n1*v3];


                N3 = @(v3) [ DN3_13*v3 - v3(end).*cos(DN3_09*v3)
                    DN3_16*v3 - v3(end).*sin(DN3_09*v3)
                    DN3_19*v3 + v3(end).*sin(DN3_09*v3)./(DN3_03*v3) - kappa*v3(end).*DN3_06*v3
                    dTN1_7p1*v1 - dTN3_9n1*v3
                    dTN3_3n1*v3 - (a+d)
                    dTN1_4p1*v1 - dTN3_6n1*v3
                    dTN3_3p1*v3 - b
                    dTN3_9p1*v3 - psib ];

            else

                N1 = @(v1) [ DN1_11*v1 - v1(end).*cos(DN1_07*v1)
                    DN1_14*v1 - v1(end).*sin(DN1_07*v1)
                    (DN1_01*v1).*(DN1_17*v1) + v1(end).*sin(DN1_07*v1) - kappa*v1(end).*(DN1_04*v1).*(DN1_01*v1)
                    dTN1_1n1*v1 - a
                    dTN1_1p1*v1 - (a+d)
                    dTN1_7n1*v1 - psia
                    dTN1_7p1*v1 - dTN3_9n1*v3
                    dTN1_4p1*v1 - dTN3_6n1*v3];


                N3 = @(v3) [ DN3_13*v3 - v3(end).*cos(DN3_09*v3)
                    DN3_16*v3 - v3(end).*sin(DN3_09*v3)
                    (DN3_03*v3).*(DN3_19*v3) + v3(end).*sin(DN3_09*v3) - kappa*v3(end).*(DN3_06*v3).*(DN3_03*v3)
                    dTN1_7p1*v1 - dTN3_9n1*v3
                    dTN3_3n1*v3 - (a+d)
                    dTN1_4p1*v1 - dTN3_6n1*v3
                    dTN3_3p1*v3 - b
                    dTN3_9p1*v3 - psib ];

            end

            % the relative norm of the residual

            res_bvp = norm(N(v),'fro')/norm(v,'fro');

            res_bvp1 = norm(N1(v1),'fro')/norm(v1,'fro');
            res_bvp3 = norm(N3(v3),'fro')/norm(v3,'fro');

            % adaptively add more Chebyshev points and resample the state
            % of the problem if the tolerance is not met.
            % Additionally, if there is too much numerical oscillation, add
            % more Chebyshev points and resample the state
            % of the problem and reset the residual.
            L21 = diffmat([2*n1 n1],0,X);
            L23 = diffmat([2*n3 n3],0,X);

            c1 = length(find(diff(sign(diff((L21*v1(2*n1 + 1:3*n1)))))));
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
                [v] = compressor_med(v1, v3);

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
                [v] = compressor_med(v1, v3);
            end

            % if the function exceeds the maximum number of iterations,
            % break with an error statement.

            if n1 > max_iter_bvp
                disp('Maximum number of Chebyshev points reached')
                res_bvp1 = res_bvp1
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
                [v] = compressor_med(v1, v3);

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
                [v] = compressor_med(v1, v3);
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

            if sqrt(kappa)*b > 10
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
                L31(1:n1-1,1:n1) = spdiags(D17*v - kappa*v(end-2).*(D04*v),0,n1-1,n1-1)*D0;
                L31((n1-1)+1:((n1-1)+(n2-1)),n1+1:n1+n2) = spdiags(D18*v - kappa*v(end-1).*(D05*v),0,n2-1,n2-1)*D00;
                L31((n1-1)+(n2-1)+1:((n1-1)+(n2-1)+(n3-1)),n1+n2+1:n1+n2+n3) = spdiags(-v(end)*sin(D09*v)./((D03*v).^2),0,n3-1,n3-1)*D000;

                L32 = spalloc((n1-1)+(n2-1)+(n3-1),n1+n2+n3,((n1-1)*n1)+((n2-1)*n2)+((n3-1)*n3));
                L32(1:n1-1,1:n1) = spdiags(-kappa*v(end-2)*(D01*v),0,n1-1,n1-1)*D0;
                L32((n1-1)+1:((n1-1)+(n2-1)),n1+1:n1+n2) = spdiags(-kappa*v(end-1)*(D02*v),0,n2-1,n2-1)*D00;
                L32((n1-1)+(n2-1)+1:((n1-1)+(n2-1)+(n3-1)),n1+n2+1:n1+n2+n3) = spdiags(-kappa*v(end))*D000;

                L33 = spalloc((n1-1)+(n2-1)+(n3-1),n1+n2+n3,((n1-1)*n1)+((n2-1)*n2)+((n3-1)*n3));
                L33(1:n1-1,1:n1) = spdiags(D01*v,0,n1-1,n1-1)*D1 + spdiags(v(end-2)*cos(D07*v),0,n1-1,n1-1)*D0;
                L33((n1-1)+1:((n1-1)+(n2-1)),n1+1:n1+n2) = spdiags(D02*v,0,n2-1,n2-1)*D10 + spdiags(v(end-1)*cos(D08*v),0,n2-1,n2-1)*D00;
                L33((n1-1)+(n2-1)+1:((n1-1)+(n2-1)+(n3-1)),n1+n2+1:n1+n2+n3) = D100 + (spdiags(v(end)*cos(D09*v)./(D03*v),0,n3-1,n3-1))*D000;

                L34 = spalloc((n1-1)+(n2-1)+(n3-1),3,((n1-1)+(n2-1)+(n3-1))*1);
                L34(1:n1-1,1:1) = sin(D07*v) - kappa*(D04*v).*(D01*v);
                L34((n1-1)+1:((n1-1)+(n2-1)),2:2) = sin(D08*v) - kappa*(D05*v).*(D02*v);
                L34((n1-1)+(n2-1)+1:((n1-1)+(n2-1)+(n3-1)),3:3) = (sin(D09*v)./(D03*v) - kappa*D06*v);

                N = @(v) [ D11*v - v(end-2).*cos(D07*v)
                    D12*v - v(end-1).*cos(D08*v)
                    D13*v - v(end).*cos(D09*v)
                    D14*v - v(end-2).*sin(D07*v)
                    D15*v - v(end-1).*sin(D08*v)
                    D16*v - v(end).*sin(D09*v)
                    (D01*v).*(D17*v) + v(end-2).*sin(D07*v) - kappa*v(end-2).*(D04*v).*(D01*v)
                    (D02*v).*(D18*v) + v(end-1).*sin(D08*v) - kappa*v(end-1).*(D05*v).*(D02*v)
                    D19*v + v(end).*sin(D09*v)./(D03*v) - kappa*v(end).*D06*v
                    dT1n1*v - a
                    dT1p1*v - dT2n1*v
                    dT7n1*v - psia
                    dT7p1*v - psic
                    dT4p1*v - dT5n1*v
                    dT2p1*v - (a + d)
                    dT8n1*v - psic
                    dT8p1*v - dT9n1*v
                    dT3n1*v - (a + d)
                    dT5p1*v - dT6n1*v
                    dT3p1*v - b
                    dT9p1*v - psib ];

                L = @(v) [ L11, L12, L13, L14
                    L21, L22, L23, L24
                    L31, L32, L33, L34
                    dT1n1
                    dT1p1 - dT2n1
                    dT7n1 
                    dT7p1
                    dT4p1 - dT5n1
                    dT2p1
                    dT8n1
                    dT8p1 - dT9n1
                    dT3n1
                    dT5p1 - dT6n1
                    dT3p1
                    dT9p1 ];
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
                L31(1:n1-1,1:n1) = spdiags(D17*v - kappa*v(end-2).*(D04*v),0,n1-1,n1-1)*D0;
                L31((n1-1)+1:((n1-1)+(n2-1)),n1+1:n1+n2) = spdiags(D18*v - kappa*v(end-1).*(D05*v),0,n2-1,n2-1)*D00;
                L31((n1-1)+(n2-1)+1:((n1-1)+(n2-1)+(n3-1)),n1+n2+1:n1+n2+n3) = spdiags(D19*v - kappa*v(end).*(D06*v),0,n3-1,n3-1)*D000;

                L32 = spalloc((n1-1)+(n2-1)+(n3-1),n1+n2+n3,((n1-1)*n1)+((n2-1)*n2)+((n3-1)*n3));
                L32(1:n1-1,1:n1) = spdiags(-kappa*v(end-2)*(D01*v),0,n1-1,n1-1)*D0;
                L32((n1-1)+1:((n1-1)+(n2-1)),n1+1:n1+n2) = spdiags(-kappa*v(end-1)*(D02*v),0,n2-1,n2-1)*D00;
                L32((n1-1)+(n2-1)+1:((n1-1)+(n2-1)+(n3-1)),n1+n2+1:n1+n2+n3) = spdiags(-kappa*v(end)*(D03*v),0,n3-1,n3-1)*D000;

                L33 = spalloc((n1-1)+(n2-1)+(n3-1),n1+n2+n3,((n1-1)*n1)+((n2-1)*n2)+((n3-1)*n3));
                L33(1:n1-1,1:n1) = spdiags(D01*v,0,n1-1,n1-1)*D1 + spdiags(v(end-2)*cos(D07*v),0,n1-1,n1-1)*D0;
                L33((n1-1)+1:((n1-1)+(n2-1)),n1+1:n1+n2) = spdiags(D02*v,0,n2-1,n2-1)*D10 + spdiags(v(end-1)*cos(D08*v),0,n2-1,n2-1)*D00;
                L33((n1-1)+(n2-1)+1:((n1-1)+(n2-1)+(n3-1)),n1+n2+1:n1+n2+n3) = spdiags(D03*v,0,n3-1,n3-1)*D100 + spdiags(v(end)*cos(D09*v),0,n3-1,n3-1)*D000;

                L34 = spalloc((n1-1)+(n2-1)+(n3-1),3,((n1-1)+(n2-1)+(n3-1))*1);
                L34(1:n1-1,1:1) = sin(D07*v) - kappa*(D04*v).*(D01*v);
                L34((n1-1)+1:((n1-1)+(n2-1)),2:2) = sin(D08*v) - kappa*(D05*v).*(D02*v);
                L34((n1-1)+(n2-1)+1:((n1-1)+(n2-1)+(n3-1)),3:3) = sin(D09*v) - kappa*(D06*v).*(D03*v);

                N = @(v) [ D11*v - v(end-2).*cos(D07*v)
                    D12*v - v(end-1).*cos(D08*v)
                    D13*v - v(end).*cos(D09*v)
                    D14*v - v(end-2).*sin(D07*v)
                    D15*v - v(end-1).*sin(D08*v)
                    D16*v - v(end).*sin(D09*v)
                    (D01*v).*(D17*v) + v(end-2).*sin(D07*v) - kappa*v(end-2).*(D04*v).*(D01*v)
                    (D02*v).*(D18*v) + v(end-1).*sin(D08*v) - kappa*v(end-1).*(D05*v).*(D02*v)
                    (D03*v).*(D19*v) + v(end).*sin(D09*v) - kappa*v(end).*(D06*v).*(D03*v)
                    dT1n1*v - a
                    dT1p1*v - dT2n1*v
                    dT7n1*v - psia
                    dT7p1*v - psic
                    dT4p1*v - dT5n1*v
                    dT2p1*v - (a + d)
                    dT8n1*v - psic
                    dT8p1*v - dT9n1*v
                    dT3n1*v - (a + d)
                    dT5p1*v - dT6n1*v
                    dT3p1*v - b
                    dT9p1*v - psib ];

                L = @(v) [ L11, L12, L13, L14
                    L21, L22, L23, L24
                    L31, L32, L33, L34
                    dT1n1
                    dT1p1 - dT2n1
                    dT7n1 
                    dT7p1
                    dT4p1 - dT5n1
                    dT2p1
                    dT8n1
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
                newt_numb = newt_numb+1;

                if kk > max_iter_newton
                    disp('Maximum number of Newton iterations reached')
                    res_newton = res_newton
                    break
                end

            end

            %% residual and other loop conditions

            [v1, v2, v3] = extractor(v, n1, n2, n3);

            if sqrt(kappa)*b <= 10

                N1 = @(v1) [ DN1_11*v1 - v1(end).*cos(DN1_07*v1)
                    DN1_14*v1 - v1(end).*sin(DN1_07*v1)
                    (DN1_01*v1).*(DN1_17*v1) + v1(end).*sin(DN1_07*v1) - kappa*v1(end).*(DN1_04*v1).*(DN1_01*v1)
                    dTN1_1n1*v1 - a
                    dTN1_1p1*v1 - dTN2_2n1*v2
                    dTN1_7n1*v1 - psia
                    dTN1_7p1*v1 - psic
                    dTN1_4p1*v1 - dTN2_5n1*v2];

                N2 = @(v2) [ DN2_12*v2 - v2(end).*cos(DN2_08*v2)
                    DN2_15*v2 - v2(end).*sin(DN2_08*v2)
                    (DN2_02*v2).*(DN2_18*v2) + v2(end).*sin(DN2_08*v2) - kappa*v2(end).*(DN2_05*v2).*(DN2_02*v2)
                    dTN2_8n1*v2 - psic
                    dTN2_2n1*v2 - dTN1_1p1*v1
                    dTN1_4p1*v1 - dTN2_5n1*v2
                    dTN2_2p1*v2 - (a + d)
                    dTN2_8p1*v2 - dTN3_9n1*v3
                    dTN2_5p1*v2 - dTN3_6n1*v3];

                N3 = @(v3) [ DN3_13*v3 - v3(end).*cos(DN3_09*v3)
                    DN3_16*v3 - v3(end).*sin(DN3_09*v3)
                    (DN3_03*v3).*(DN3_19*v3) + v3(end).*sin(DN3_09*v3) - kappa*v3(end).*(DN3_06*v3).*(DN3_03*v3)
                    dTN2_8p1*v2 - dTN3_9n1*v3
                    dTN3_3n1*v3 - (a + d)
                    dTN2_5p1*v2 - dTN3_6n1*v3
                    dTN3_3p1*v3 - b
                    dTN3_9p1*v3 - psib ];

            else
                 N1 = @(v1) [ DN1_11*v1 - v1(end).*cos(DN1_07*v1)
                    DN1_14*v1 - v1(end).*sin(DN1_07*v1)
                    (DN1_01*v1).*(DN1_17*v1) + v1(end).*sin(DN1_07*v1) - kappa*v1(end).*(DN1_04*v1).*(DN1_01*v1)
                    dTN1_1n1*v1 - a
                    dTN1_1p1*v1 - dTN2_2n1*v2
                    dTN1_7n1*v1 - psia
                    dTN1_7p1*v1 - psic
                    dTN1_4p1*v1 - dTN2_5n1*v2];

                N2 = @(v2) [ DN2_12*v2 - v2(end).*cos(DN2_08*v2)
                    DN2_15*v2 - v2(end).*sin(DN2_08*v2)
                    (DN2_02*v2).*(DN2_18*v2) + v2(end).*sin(DN2_08*v2) - kappa*v2(end).*(DN2_05*v2).*(DN2_02*v2)
                    dTN2_8n1*v2 - psic
                    dTN2_2n1*v2 - dTN1_1p1*v1
                    dTN1_4p1*v1 - dTN2_5n1*v2
                    dTN2_2p1*v2 - (a + d)
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

%% The function for exctracting v1, v2, v3

    function [v1, v3] = extractor_med(v, n1, n3)

        v1 = zeros(3*n1+1,1);
        v1(1:n1) = v(1:n1);
        v1(n1+1:2*n1) = v(n1+n3+1:2*n1+n3);
        v1(2*n1+1:3*n1) = v(2*n1+2*n3+1:3*n1+2*n3);
        v1(end) = v(end-1);

        v3 = zeros(3*n3+1,1);
        v3(1:n3) = v(n1+1:n1+n3);
        v3(n3+1:2*n3) = v(2*n1+n3+1:2*n1+2*n3);
        v3(2*n3+1:3*n3) = v(3*n1+2*n3+1:3*n1+3*n3);
        v3(end) = v(end);

    end

%% The function for compressing v1, v2, v3

    function [v] = compressor_med(v1, v3)

        cn1 = (length(v1) - 1)/3;
        cn3 = (length(v3) - 1)/3;
        v = zeros(3*(cn1+cn3)+2,1);
        v(1:cn1) = v1(1:cn1);
        v(cn1+1:cn1+cn3) = v3(1:cn3);
        v(cn1+cn3+1:2*cn1+cn3) = v1(cn1+1:2*cn1);
        v(2*cn1+cn3+1:2*cn1+2*cn3) = v3(cn3+1:2*cn3);
        v(2*cn1+2*cn3+1:3*cn1+2*cn3) = v1(2*cn1+1:3*cn1);
        v(3*cn1+2*cn3+1:3*cn1+3*cn3) = v3(2*cn3+1:3*cn3);
        v(end-1) = v1(end);
        v(end) = v3(end);




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

if ((abs(psib_actual) > pi/2)||(abs(psia_actual) > pi/2))
    %% Assigning the output into variables that have a physical meaning
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
    plot([a;a],ylim,'-.k')
    plot([b;b],ylim,'-.k')
    hold off
else
    %% Assigning the output into variables that have a physical meaning

    R1 = v(1:n1);
    R3 = v(n1+1:n1+n3);
    U1 = v(n1+n3+1:2*n1+n3);
    U3 = v(2*n1+n3+1:2*n1+2*n3);
    Psi1 = v(2*n1+2*n3+1:3*n1+2*n3);
    Psi3 = v(3*n1+2*n3+1:3*n1+3*n3);
    ell01 = v(end-1);
    ell03 = v(end);

    %% Plotting
    % delete or comment this block if the function is to be used inside some
    % other algortihm.
    figure(3)
    plot(chebfun(R1),chebfun(U1),'.-k')
    hold on
    plot(chebfun(R3),chebfun(U3),'.-k')
    axis equal
    plot([a;a],ylim,'-.k')
    plot([b;b],ylim,'-.k')
    hold off
end


%% Tracking n_v and n_N
    n_v = length(v);
    n_N = newt_numb;
end