function sol = leapfrog(x0, v0, tspan)
	n = length(tspan);
	sol = zeros(1, n);	
	sol(1) = x0;

	h = tspan(2) - tspan(1);
	x = x0;
	vhalf = 0;
	v = v0;
	
	for i=2:n
		vhalf = v + 0.5 * h * -x;
		x = x + h * vhalf;
		v = vhalf + 0.5 * h * -x;
		sol(i) = x;
	end	
end
