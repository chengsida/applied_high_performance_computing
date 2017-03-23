function RTstability()
    clear all; close all; clc;
    [X, Y]          = meshgrid(-4:0.1:4, -4:0.1:4);
    Z               = X + i*Y;
    sigma           = abs(1 + Z + (Z.^2)/2 + (Z.^3)/6 + (Z.^4)/24);
    figure('WindowStyle', 'docked');
    contourf(X, Y, sigma, [1 1],'-k');
    hold on;
    axis('equal', [-4 4 -4 4]);
    grid on;
    xlabel('\lambda_{Re}\Delta t');
    ylabel('\lambda_{Im}\Delta t');