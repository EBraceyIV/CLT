% EAS4240
% Final Exam, Q4
% Edgar Bracey IV
clc; clear;

% % Mechanical properties
E1 = 155; E2 = 12.1; G12 = 4.4; % GPa
E1 = E1*1E9; E2 = E2*1E9; G12 = G12*1E9;
nu12 = 0.248; nu21 = nu12/E1 * E2;
t = 0.15; % mm
t = t/1000;
layers = [45; 90; -45; -45; 90; 45];
N = length(layers);

% Initialize S & Q matrices
S = zeros(3,3,N);
Q = zeros(3,3,N);
Qbar = zeros(3,3,N);
Sbar = zeros(3,3,N);

% Build S & Q matrices
for layer = 1:N
    S(1,1,layer) = 1/E1;
    S(1,2,layer) = -nu12/E1;
    S(2,1,layer) = -nu21/E2;
    S(2,2,layer) = 1/E2;
    S(3,3,layer) = 1/G12;
    Q(:,:,layer) = inv(S(:,:,layer));
end

% Build Sbar and Qbar matrices
for layer = 1:N
    theta = layers(layer);
    c = cosd(theta); s = sind(theta);
    Qbar(1,1,layer) = Q(1,1,layer)*c^4 + Q(2,2,layer)*s^4 + 2*(Q(1,2,layer) + 2*Q(3,3,layer))*s^2*c^2;
    Qbar(1,2,layer) = (Q(1,1,layer) + Q(2,2,layer) - 4*Q(3,3,layer))*s^2*c^2 + Q(1,2,layer)*(c^4 + s^4);
    Qbar(2,2,layer) = Q(1,1,layer)*s^4 + Q(2,2,layer)*c^4 + 2*(Q(1,2,layer) + 2*Q(3,3,layer))*s^2*c^2;
    Qbar(1,3,layer) = (Q(1,1,layer) - Q(1,2,layer) - 2*Q(3,3,layer))*c^3*s - (Q(2,2,layer) - Q(1,2,layer) - 2*Q(3,3,layer))*c*s^3;
    Qbar(2,3,layer) = (Q(1,1,layer) - Q(1,2,layer) - 2*Q(3,3,layer))*c*s^3 - (Q(2,2,layer) - Q(1,2,layer) - 2*Q(3,3,layer))*c^3*s;
    Qbar(3,3,layer) = (Q(1,1,layer) + Q(2,2,layer) - 2*Q(1,2,layer) - 2*Q(3,3,layer))*s^2*c^2 + Q(3,3,layer)*(s^4 + c^4);
    Qbar(2,1,layer) = Qbar(1,2,layer);
    Qbar(3,1,layer) = Qbar(1,3,layer);
    Qbar(3,2,layer) = Qbar(2,3,layer);
    Sbar(:,:,layer) = inv(Qbar(:,:,layer));
end

% Calculate the range of z's for the layers
z = zeros(1, N + 1);
for ii = 1:length(z)
    z(ii) = (N*t)/2 - t*(ii - 1);
end

% Calculate the ABD matrices
A = zeros(3,3); B = zeros(3,3); D = zeros(3,3);
for ii = 1:3
    for jj = 1:3
        for kk = 1:N
            A(ii,jj) = A(ii,jj) + Qbar(ii,jj,kk) * (z(kk) - z(kk + 1));
            B(ii,jj) = B(ii,jj) + Qbar(ii,jj,kk) * ((z(kk))^2 - (z(kk + 1))^2);
            D(ii,jj) = D(ii,jj) + Qbar(ii,jj,kk) * ((z(kk))^3 - (z(kk + 1))^3);
        end
    end
end
B = -B/2; D = D/3;

% Build ABD matrix
ABD = [A B; B D];

% Invert the ABD matrices
ABDi = inv(ABD);

% Failure properties [Pa]
SLp = 1500E6; SLn = -1250E6;
STp = 50E6; STn = -200E6;
SLT = 100E6;

% Tsai-Wu F values
F11 = 1/(SLp*abs(SLn));
F22 = 1/(STp*abs(STn));
F66 = 1/SLT^2;
F1 = 1/SLp - 1/abs(SLn);
F2 = 1/STp - 1/abs(STn);

% This range was found by manually testing every 10000 N/m until failure
for load = 90000:1:100000
    NM = [load; 0; 0; 0; 0; 0];
    epKap = ABDi*NM;
    
    % Initalize transformation and global & local stress/strain matrices
    % Value at top of layer is column 1 (left), bottom is column 2 (right)
    strain12 = zeros(3,2,N);
    strainxy = zeros(3,2,N);
    stress12 = zeros(3,2,N);
    stressxy = zeros(3,2,N);
    
    % Transformation matrices
    Teps = zeros(3,3,N);
    Tsig = zeros(3,3,N);
    
    for ii = 0:N-1
        theta = layers(ii+1);
        c = cosd(theta); s = sind(theta);
        
        Teps = [c^2 s^2 -c*s; s^2 c^2 -c*s; -2*c*s 2*c*s (c^2 - s^2)];
        Tsig = [c^2 s^2 2*c*s; s^2 c^2 -2*c*s; -c*s c*s (c^2 - s^2)];
        
        % Global coordinates
        strainxy(:,1,ii+1) = epKap(1:3) + z(length(z) - ii)*epKap(4:6);% - deltaT*[alphax; alphay; alphaxy];%alphaxy(:,1,ii+1);
        strainxy(:,2,ii+1) = epKap(1:3) + z(length(z) - ii - 1)*epKap(4:6);% - deltaT*[alphax; alphay; alphaxy];%alphaxy(:,1,ii+1);
        stressxy(:,1,ii+1) = Qbar(:,:,ii+1)*strainxy(:,1,ii+1);
        stressxy(:,2,ii+1) = Qbar(:,:,ii+1)*strainxy(:,2,ii+1);
        
        % Local coordinates
        stress12(:,1,ii+1) = Tsig*stressxy(:,1,ii+1);
        stress12(:,2,ii+1) = Tsig*stressxy(:,2,ii+1);
    end

    for ii = 1:N
        % Coefficients of the polynomial where the root, lambda, 
        % is the factor of safety
        TW(1,:,ii) = [(F11*stress12(1,1,ii)^2 + F22*stress12(2,1,ii)^2 + F66*stress12(3,1,ii)^2 - sqrt(F11*F22)*stress12(1,1,ii)*stress12(2,1,ii)) (F1*stress12(1,1,ii) + F2*stress12(2,1,ii)) -1];
        TW(2,:,ii) = [(F11*stress12(1,2,ii)^2 + F22*stress12(2,2,ii)^2 + F66*stress12(3,2,ii)^2 - sqrt(F11*F22)*stress12(1,2,ii)*stress12(2,2,ii)) (F1*stress12(1,2,ii) + F2*stress12(2,2,ii)) -1];
        
        lambda(1,:,ii) = roots(TW(1,:,ii));
        lambda(2,:,ii) = roots(TW(2,:,ii));
    end
    
    % This finds the minimum factor of safety in the positive range
    % between the bottom and top of each layer
    SF = min(max(max(lambda)));
    
    % When the factor of safety is below 1, the laminate experiences the
    % first ply failure at the given load
    if SF < 1
        fprintf('First-ply-failure occurs at Nx = %d N/m\n',load)
        break;
    end
end