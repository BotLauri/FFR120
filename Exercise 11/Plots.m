function Plots(lattice, SIR, beta, gamma, t)
    subplot(1, 2, 1)
    PlotLattice(lattice)
    title(strcat('t = ', int2str(t)))
    subplot(1, 2, 2)
    hold on
    plot(0:t, SIR(:, 1), 'Color', 'b')
    plot(0:t, SIR(:, 2), 'Color', '#FFA500')
    plot(0:t, SIR(:, 3), 'Color', 'g')
    title(strcat(' \beta = ', num2str(beta, 2)))
    title(strcat(' \gamma = ', num2str(gamma, 2)))
    pbaspect([1 1 1])
    legend('S', 'I', 'R')
    hold off
end

