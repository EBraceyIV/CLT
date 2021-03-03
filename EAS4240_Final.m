% EAS4240
% Final Exam
% Edgar Bracey IV
clc; clear;

% % Mechanical properties
E1 = 155; E2 = 12.1; G12 = 4.4; % GPa
E1 = E1*1E9; E2 = E2*1E9; G12 = G12*1E9;
nu12 = 0.248; nu21 = nu12/E1 * E2;
t = 0.15; % mm
t = t/1000;

% Build laminate
N = input('Enter number of layers: ');
layers = zeros(N,1);
for ii = 1:N
    fprintf('Enter fiber diretion of layer %g: ',ii)
    layers(ii) = input('');
end

% Thermal properties
alpha1 = -0.018E-6; alpha2 = 24.3E-6; % 10^-6 m/m/C
alpha12 = [alpha1; alpha2; 0];
thermBool = input('Apply thermal effects? (y/n) ','s');
if (thermBool == 'y')
    deltaT = input('Enter change in temperature (C): ');
else
    deltaT = 0; %C
end

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
Ai = ABDi(1:3,1:3);
Bi = ABDi(1:3,4:6);
Di = ABDi(4:6,4:6);

%Effective constants for the laminate
Ex = 1/(t*ABDi(1,1));
Ey = 1/(t*ABDi(2,2));
Gxy = 1/(t*ABDi(3,3));
nuxy = -ABDi(1,2)/ABDi(1,1);

% Decide to apply stress or strain
dims = {'x' 'y' 'xy'};
NM = zeros(6,1);
epKap = zeros(6,1);
appChoice = 2;
while appChoice ~= 1 || appChoice ~= 2
        fprintf('Choose 1 to apply stress, or 2 to apply strain: ')
        appChoice = input('');
    if appChoice == 1 % Apply force & moment resultants, return mid-plane strain and curvature
                for NN = 1:3
                    fprintf('Enter force resultant N%s [N/m]: ',dims{NN})
                    NM(NN) = input('');
                end
                for MM = 4:6
                    fprintf('Enter moment resultant M%s [Nm/m]: ',dims{MM-3})
                    NM(MM) = input('');
                end
        epKap = ABDi*NM;
        break
    elseif appChoice ==2 % Apply mid-plane strain and curvature, return force & moment resultants
        for eps = 1:3
            fprintf('Enter mid-plane strain epsilon0%s [microstrain]: ',dims{eps})
            epKap(eps) = input('');
            epKap(eps) = epKap(eps)/1E6;
        end
        for kap = 4:6
            fprintf('Enter curvatuve k%s [m^-1]: ',dims{kap-3})
            epKap(kap) = input('');
        end
        NM = ABD*epKap;
        break
    end
end

% Initalize transformation and global & local stress/strain matrices 
% Value at top of layer is column 1 (left), bottom is column 2 (right)
strain12 = zeros(3,2,N);
strainxy = zeros(3,2,N);
stress12 = zeros(3,2,N);
stressxy = zeros(3,2,N);

% Transformation matrices
Teps = zeros(3,3,N);
Tsig = zeros(3,3,N);

% Thermal force resultants initialized
NxT = 0; NyT = 0; NxyT = 0;

for ii = 0:N-1
    theta = layers(ii+1);
    c = cosd(theta); s = sind(theta);
    
    % Transformation matrices
    Teps = [c^2 s^2 -c*s; s^2 c^2 -c*s; -2*c*s 2*c*s (c^2 - s^2)];
    Tsig = [c^2 s^2 2*c*s; s^2 c^2 -2*c*s; -c*s c*s (c^2 - s^2)];
    
    % Structural thermal coefficients
    alphax = alpha1*c^2 + alpha2*s^2;
    alphay = alpha1*s^2 + alpha2*c^2;
    alphaxy = 2*(alpha1 - alpha2)*c*s;
    
    % Thermal Force Resultants
    NxT = NxT + Qbar(1,:,ii+1) * deltaT*[alphax; alphay; alphaxy] * (z(ii + 1) - z(ii + 2));
    NyT = NyT + Qbar(2,:,ii+1) * deltaT*[alphax; alphay; alphaxy] * (z(ii + 1) - z(ii + 2));
    NxyT = NxyT + Qbar(3,:,ii+1) * deltaT*[alphax; alphay; alphaxy] * (z(ii + 1) - z(ii + 2));
    
    % Global coordinates
    strainxy(:,1,ii+1) = epKap(1:3) + z(length(z) - ii)*epKap(4:6);% - deltaT*[alphax; alphay; alphaxy];
    strainxy(:,2,ii+1) = epKap(1:3) + z(length(z) - ii - 1)*epKap(4:6);% - deltaT*[alphax; alphay; alphaxy];
    stressxy(:,1,ii+1) = Qbar(:,:,ii+1)*(strainxy(:,1,ii+1) - deltaT*[alphax; alphay; alphaxy]);
    stressxy(:,2,ii+1) = Qbar(:,:,ii+1)*(strainxy(:,2,ii+1) - deltaT*[alphax; alphay; alphaxy]);
    
    % Local coordinates
    strain12(:,1,ii+1) = Teps*strainxy(:,1,ii+1);
    strain12(:,2,ii+1) = Teps*strainxy(:,2,ii+1);
    stress12(:,1,ii+1) = Tsig*stressxy(:,1,ii+1);
    stress12(:,2,ii+1) = Tsig*stressxy(:,2,ii+1);
end

% Laminate Coeffs. of Thermal Expansion
aXhat = ABDi(1,1:3)*[NxT; NyT; NxyT]/deltaT;
aYhat = ABDi(2,1:3)*[NxT; NyT; NxyT]/deltaT;
aXYhat = ABDi(3,1:3)*[NxT; NyT; NxyT]/deltaT;

% Make some cool graphs?
doPlot = input("Plot? (y/n) ",'s');
if (doPlot == 'y')
    % Plot the structural strains
    plotstrainx = [];
    plotstrainy = [];
    plotstrainxy = [];
    for ii = 1:N
        addLayerx = linspace(strainxy(1,1,ii),strainxy(1,2,ii));
        plotstrainx = [plotstrainx addLayerx];
        addLayery = linspace(strainxy(2,1,ii),strainxy(2,2,ii));
        plotstrainy = [plotstrainy addLayery];
        addLayerxy = linspace(strainxy(3,1,ii),strainxy(3,2,ii));
        plotstrainxy = [plotstrainxy addLayerxy];
    end
    
    figure
    plot(plotstrainx*1E6,linspace(z(end),z(1),length(plotstrainx))*1000,'-','LineWidth',2)
    hold on
    plot(plotstrainy*1E6,linspace(z(end),z(1),length(plotstrainy))*1000,'-','LineWidth',2)
    hold on
    plot(plotstrainxy*1E6,linspace(z(end),z(1),length(plotstrainxy))*1000,'-','LineWidth',2)
    xlabel('Structural Strains [\mu\epsilon]')
    ylabel('z [mm]')
    legend('\epsilon_x','\epsilon_y','\gamma_{xy}')
    set(gca, 'YDir','reverse'); grid on
    
    % Plot the principal material strains
    plotstrain1 = [];
    plotstrain2 = [];
    plotstrain12 = [];
    for ii = 1:N
        addLayer1 = linspace(strain12(1,1,ii),strain12(1,2,ii));
        plotstrain1 = [plotstrain1 addLayer1];
        addLayer2 = linspace(strain12(2,1,ii),strain12(2,2,ii));
        plotstrain2 = [plotstrain2 addLayer2];
        addLayer12 = linspace(strain12(3,1,ii),strain12(3,2,ii));
        plotstrain12 = [plotstrain12 addLayer12];
    end
    
    figure
    plot(plotstrain1*1E6,linspace(z(end),z(1),length(plotstrain1))*1000,'-','LineWidth',2)
    hold on
    plot(plotstrain2*1E6,linspace(z(end),z(1),length(plotstrain2))*1000,'-','LineWidth',2)
    hold on
    plot(plotstrain12*1E6,linspace(z(end),z(1),length(plotstrain12))*1000,'-','LineWidth',2)
    xlabel('Principal Material Strains [\mu\epsilon]')
    ylabel('z [mm]')
    legend('\epsilon_1','\epsilon_2','\gamma_{12}')
    set(gca, 'YDir','reverse'); grid on
    
    % Plot the structural stresses
    plotstressx = [];
    plotstressy = [];
    plotstressxy = [];
    for ii = 1:N
        addLayerx = linspace(stressxy(1,1,ii),stressxy(1,2,ii));
        plotstressx = [plotstressx addLayerx];
        addLayery = linspace(stressxy(2,1,ii),stressxy(2,2,ii));
        plotstressy = [plotstressy addLayery];
        addLayerxy = linspace(stressxy(3,1,ii),stressxy(3,2,ii));
        plotstressxy = [plotstressxy addLayerxy];
    end
    
    figure
    plot(plotstressx/1E6,linspace(z(end),z(1),length(plotstressx))*1000,'-','LineWidth',2)
    hold on
    plot(plotstressy/1E6,linspace(z(end),z(1),length(plotstressy))*1000,'-','LineWidth',2)
    hold on
    plot(plotstressxy/1E6,linspace(z(end),z(1),length(plotstressy))*1000,'-','LineWidth',2)
    xlabel('Structual Stresses [MPa]')
    ylabel('z [mm]')
    legend('\sigma_x','\sigma_y','\tau_{xy}')
    set(gca, 'YDir','reverse'); grid on
    
    % Plot the principal material stresses
    plotstress1 = [];
    plotstress2 = [];
    plotstress12 = [];
    for ii = 1:N
        addLayer1 = linspace(stress12(1,1,ii),stress12(1,2,ii));
        plotstress1 = [plotstress1 addLayer1];
        addLayer2 = linspace(stress12(2,1,ii),stress12(2,2,ii));
        plotstress2 = [plotstress2 addLayer2];
        addLayer12 = linspace(stress12(3,1,ii),stress12(3,2,ii));
        plotstress12 = [plotstress12 addLayer12];
    end
    
    figure
    plot(plotstress1/1E6,linspace(z(end),z(1),length(plotstress1))*1000,'-','LineWidth',2)
    hold on
    plot(plotstress2/1E6,linspace(z(end),z(1),length(plotstress2))*1000,'-','LineWidth',2)
    hold on
    plot(plotstress12/1E6,linspace(z(end),z(1),length(plotstress12))*1000,'-','LineWidth',2)
    xlabel('Principal Material Stresses [MPa]')
    ylabel('z [mm]')
    legend('\sigma_1','\sigma_2','\tau_{12}')
    set(gca, 'YDir','reverse'); grid on
end