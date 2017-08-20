%% 1D DIC code by Thorsten Becker
close all; clear all; clc;
 
%% Generate a random function F with noise (i.e. 1D image) - sin(Ax+B) where B = displacment and A = stretch.
X = 0:1:20;
F = sin(X/5) +rand(size(X)).*0.3;      % zero load image (with noise).
u = [0.2; 1.4];                         % known displacment and stretch - this is what we are trying to slove for.
G = sin(u(2).*X/5+u(1)) +rand(size(X)).*0.3;  % impose known displacment and strech on F (with more noise).
 
%% Choose subset
subsize = 15;                           % subset size
offset = 5;                             % subset position
f = F(offset:offset+subsize-1);         % subset intensity values
x = X(offset:offset+subsize-1);         % subset co-ordinates
x0 = subsize/2 + offset;                % subset centre
 
%% Initialise solver - 1st order                   
dfdx = [gradient(f)', gradient(f)'.*(x-x0)'];
H = (dfdx'*dfdx);                       % G'G - pg.86 (Hesian matrix - only need to compute this once when using method on pg.100)
d = [0;0];                              % initial guess 
 
for i = 1:20;                           % solve to n iterations
    g = interp1(X,G,(x + (d(1,i)+d(2,i)*(x-x0))),'spline'); %compute values to sub-pixel accuracy, use 'nearest' for no subpixel accuracy
    q = dfdx'*(f-g)';                   % compute G'g - pg.86
    dd = pinv(H) * q;                   % calulate new displacement guess based on the difference between f and g.
    d(:,i+1) = d(:,i) + dd;             % update displacment guess.
end
 
%% plot some cool stuff
figure;
subplot(1,2,1), plot(x + d(1,end) + d(2,end)*(x-x0),f,'r-','linewidth',2); hold on
subplot(1,2,1), plot(x,f,'b-','linewidth',2);
subplot(1,2,1),plot(X,G,'r--o'); 
subplot(1,2,1),plot(X,F,'b--o'); legend('deformed subset','original subset','G','F'); hold off
subplot(1,2,2),plot(d(1,:),'-ro'); hold on;
subplot(1,2,2),plot(d(2,:),'-bo'); legend('convergene of u','convergence of du/dx'); hold off