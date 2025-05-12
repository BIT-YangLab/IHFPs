
rng('default');  % For reproducibility
x = randn(100, 1);
y = 2 * randn(100, 1) + 1;  


figure;
t = tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');


nexttile(3, [1, 1]);  
scatter(x, y, 'filled');
xlabel('X-axis');
ylabel('Y-axis');
grid on;

nexttile(1, [1, 1]);  
histogram(x, 10, 'Normalization', 'pdf'); 

xlabel('X-axis');
ylabel('Density');
grid on;


nexttile(4, [1, 1]);  
histogram(y, 10, 'Normalization', 'pdf', 'Orientation', 'horizontal');

xlabel('Density');
ylabel('Y-axis');
grid on;
