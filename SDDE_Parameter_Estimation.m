% This experiment aims to determine the drift/volatility operators for the
% one-dimensional SDDE
% dX(t) = (mu0*X(t)+mu1*X(t-tau)) dt  +  sigma dW(t),   0<t<T,    (*)
% X(t) = phi(t), t<=0,
% for given delay tau and history phi. In this experiment, we assume we
% know the value of the delay tau, and we know the SDDE takes the form (*), with
% constants mu0, mu1, and sigma. However, we do not know the values of
% these constants (mu0, mu1, and sigma), but we do have a selection of
% "number" simulated sample paths.

rng('default')

% Equation parameters:
T = 10; % Terminal time.
mu0 = -1; % First drift constant.
mu1 = 0.8; % Second drift constant.
sigma = 0.2; % Volatility constant.
tau = 2; % The time delay. Keep this an integer, for this script.
phi = @(t) -t.^3;

% Simulation parameters:
dt = 2^-7; % Time increment.
number = 100; % Number of trials.

% Simulate X with a simple Euler--Maruyama method:
t = dt:dt:T; N = T/dt; % Set the time axis.
W = cumsum(sqrt(dt)*randn(number,N),2); % Wiener paths.
Xpre = phi(-tau:dt:-dt); X0 = phi(0); % Define process history.
X = zeros(number,length(t))/2; % Initialise the vectorised solution.
Xd = zeros(number,length(t)); % Initialise the delayed solution Xd(t) = X(t-1).
% EM numerical solution:
for n = 1:N-1
    if n == 1
        X(:,n) = X0*ones(number,1);
    end
    if n <= round(tau/dt)
        Xd(:,n) = phi(t(n)-tau);
    else
        Xd(:,n) = X(:,n - round(tau/dt));
    end
    X(:,n+1) = X(:,n) + (mu0*X(:,n)+mu1*Xd(:,n))*dt + sigma*(W(:,n+1)-W(:,n));
end
Xd(:,end) = X(:,end - round(tau/dt)); % Not needed but I have this for completeness.

% Use the process increments to estimate sigma (with quadratic variation):
dX = X(:,2:end)-X(:,1:end-1);
sigma_estimate = sqrt(mean(dX(:).^2)/dt);

% Use regression (least squares) to estimate drift coefficients
Xdata  = X(:,1:end-1);   Xdata = Xdata(:);
Xddata = Xd(:,1:end-1);  Xddata = Xddata(:);
driftdata = dX(:) / dt;
A = [Xdata, Xddata];
regression = (A'*A) \ (A'*driftdata);
mu0_estimate = regression(1);
mu1_estimate = regression(2);

% Simulate the deterministic solution of dx(t)/dt = mu0*x(t)+mu1*X(t-tau):
% This could also be done with the analytic formula, in more efficient way,
% but I'm choosing to use the Euler method.
x = ones(1,length(t)); xd = ones(1,length(t));
x(1) = X0 + (mu0*X0+mu1*phi(-tau))*dt;
for n = 1:N-1
    % Note x(0)=1 already.
    td = t(n) - tau;
    if n <= round(tau/dt)
        xd(n) = phi(td);
    else
        xd(n) = x(n-round(tau/dt));
    end
    x(n+1) = x(n) + (mu0*x(n)+mu1*xd(n))*dt;
end

% Estimated deterministic solution of dy(t)/dt = mu0_estimate*y(t)+mu1_estimate*y(t-tau):
y = ones(1,length(t));
yd = ones(1,length(t));
y(1) = X0 + (mu0_estimate*X0+mu1_estimate*phi(-tau))*dt;
for n = 1:N-1
    td = t(n)-tau;
    if n <= round(tau/dt)
        yd(n) = phi(td);
    else
        yd(n) = y(n-round(tau/dt));
    end
    y(n+1) = y(n) + (mu0_estimate*y(n)+mu1_estimate*yd(n))*dt;
end

figure(1)
title(['Estimated SDDE Parameters from ',num2str(number),' Trials'],'Interpreter','latex','FontSize',16)
xlim([-1,T])
hold on
L1=plot(nan,nan,'-','Color',[0.25 0.40 0.85 0.35],'linewidth',1);
L2=plot(nan,nan,'-','Color',[0.85 0.10 0.10],'linewidth',2);
L3=plot(nan,nan,'--','Color',[0.90 0.45 0.0],'linewidth',2);
L4=plot(nan,nan,'.-','Color',[0.45 0.90 0.0],'linewidth',2);
legend([L1,L2,L3,L4],{'$X_i(t)$, $i=1,\ldots,100$','$\frac{\mathrm{d}x(t)}{\mathrm{d}t}=\mu_0x(t)+\mu_1x(t-\tau)$    ','$\frac{\mathrm{d}y(t)}{\mathrm{d}t}=\hat{\mu}_0y(t)+\hat{\mu}_1y(t-\tau)$','$y\pm\hat{\sigma}$'},'Interpreter','latex','FontSize',16)
for i = 1:min([100,number])
    plot([-tau:dt:-dt,0,t],[Xpre,phi(0),X(i,:)],'Color','[0.25 0.40 0.85 0.35]','HandleVisibility','off')
end
plot([0,t],[phi(0),x],'-','Color',[0.85 0.10 0.10],'linewidth',2,'HandleVisibility','off');
plot([0,t],[phi(0),y],'--','Color',[0.90 0.45 0.0],'linewidth',2,'HandleVisibility','off');
plot([0,t],[phi(0),y+sigma_estimate],'.-','Color',[0.45 0.90 0.0],'linewidth',2,'HandleVisibility','off');
plot([0,t],[phi(0),y-sigma_estimate],'.-','Color',[0.45 0.90 0.0],'linewidth',2,'HandleVisibility','off');
xlabel('$t$','Interpreter','latex','FontSize',16)
ylabel('$X(t)$','Interpreter','latex','FontSize',16)
set(gca,'color','[.9,.9,.9]');
set(gcf,'color','[1,1,1]');
set(gcf, 'InvertHardcopy', 'off')
set(gcf,'position',[200,200,800,800]-200)
axis([-tau,T,-0.2,3.2])

% Annotations in the figure:
param_text = sprintf(['True $\\mu_0$: %.4f\n' 'Learnt $\\hat{\\mu}_0$: %.4f\n' 'True $\\mu_1$: %.4f\n' 'Learnt $\\hat{\\mu}_1$: %.4f\n' 'True $\\sigma$: %.4f\n' 'Learnt $\\hat{\\sigma}$: %.4f'], ...
    mu0, mu0_estimate, mu1, mu1_estimate, sigma, sigma_estimate);
xlims = xlim;
ylims = ylim;
x_margin = 0.03 * (xlims(2) - xlims(1));   % 3% along horizontal margin.
y_margin = 0.4 * (ylims(2) - ylims(1));   % 40% vertical margin along.
text(xlims(2) - x_margin,ylims(1) + y_margin, param_text,'Interpreter','latex', 'FontSize',14,'HorizontalAlignment','right','VerticalAlignment','bottom','BackgroundColor','white','Margin',6,'EdgeColor','black');

% Save the figure:
if number == 100
    saveas(gcf,'SDDE_Parameter_Estimation_100.eps','epsc')
end

% Displaying the results in Matlab:
disp(' [ True Values | Estimated Values ] ')
fprintf('mu0:   %6.3f   %8.3f\n', mu0, mu0_estimate);
fprintf('mu1:   %6.3f   %8.3f\n', mu1, mu1_estimate);
fprintf('sigma: %6.3f   %8.3f\n', sigma, sigma_estimate);