function connect_dots(lambda, x, options)

options.null = 0;
col = getoptions(options, 'col', 'k');
lw = getoptions(options, 'lw', 1);
tol = getoptions(options, 'tol', 1e-1);

hold on;
for i=2:length(lambda)
    u = x{i}; v = x{i-1};
    lu = lambda(i); lv = lambda(i-1);
    for j=1:length(u)
        % find clostest neighbors
        d = abs( u(j)-v );
        I = find(d<tol);
        for k=1:length(I)
            plot([lv lu], [v(I(k)) u(j)], col, 'LineWidth', lw);
        end
    end
end


