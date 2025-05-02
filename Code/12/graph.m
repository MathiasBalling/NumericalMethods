function solve_bvp
    % Define the interval
    solinit = bvpinit(linspace(0, 2, 10), [0 0]); % Initial guess
    
    % Solve the BVP
    sol = bvp4c(@odefun, @bcfun, solinit);
    
    % Plot the solution
    x = linspace(0, 2, 1000);
    y = deval(sol, x);
    
    plot(x, y(1, :), 'b-', 'LineWidth', 2)
    xlabel('x')
    ylabel('y(x)')
    title('Solution of the BVP')
    grid on

    % Nested function: system of ODEs
    function dydx = odefun(x, y)
        dydx = [y(2); 
                2*x + sin(y(2)) - cos(y(1))];
    end

    % Nested function: boundary conditions
    function res = bcfun(ya, yb)
        res = [ya(1) - 0;  % y(0) = 0
               yb(1) - 1]; % y(2) = 1
    end
end

solve_bvp()