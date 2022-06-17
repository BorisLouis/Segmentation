function F=gaussian(x,data)

F = x(3) .* exp(-((data-x(2))./(sqrt(2).*x(1))).^2)+x(4);