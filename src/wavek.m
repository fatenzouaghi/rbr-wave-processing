function wavek=wavek(f,dep,g)
% inverts the linear dispersion relation sig^2=g*k*tanh(k*dep) to get 
% k from T=2*pi/sig
eps=0.0001;	
sig=2.*pi.*f;
Y=dep.*sig.^2./g ;   % a is the squared adimensional frequency
X=sqrt(Y);
I=1;
F=1.;
while abs(max(F)) > eps
    H=tanh(X);
    F=Y-X.*H;
    FD=-H-X./cosh(X).^2;
    X=X-F./FD;
end
wavek=X./dep;
end
