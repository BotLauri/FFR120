%% SIR model.
tic
L = 5; % Length of the lattice.
N = 50; % Number of agents. 
T = 20; % Number of timesteps. 
d = 0.8;
gamma = 0.01;

% Initialization of the grid. 
lattice = cell(L);
for n = 1:N
    xCoord = randi(L);
    yCoord = randi(L);
    lattice{xCoord, yCoord}{end + 1} = 1;
end

%imagesc(lattice)
for t = 1:T
    lattice = Diffusion(lattice, d);
    %lattice = Infection(lattice, beta);
    %lattice = Recovery(lattice, gamma);
    %imagesc(lattice)
    pause(0.2)
end

toc
