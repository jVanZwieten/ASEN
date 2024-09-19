clear

Ehat = [-.8624; -.3624; -.3535];
Nhat = [-sqrt(.5); -sqrt(.5); 0];
norm(Nhat)
omega = acos(dot(Ehat, Nhat))
omegaD = rad2deg(omega)