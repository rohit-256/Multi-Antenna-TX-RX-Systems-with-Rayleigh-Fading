function [X_0,X_k] = Jake(v,fc,k,no_samples,L)
    %v m/s
    % fc Mhz
    %k time instant
    %Ts signalling interval
    %L number of samples
    c=300000000;                    %speed of light     
    fd = fc*v/c;                    %Doppler Frequency
    wd = 2*pi*fd;
    Ts = 1/10*(1/(2*fd));           %sampling time = 0.1 * coherence time
    N = 100;                        % Number of signals
    M= N/4; 


    psi = 2*pi*rand(no_samples*L,M) - pi;
    phi = 2*pi*rand(no_samples*L,M) - pi;
    theta = 2*pi*rand(no_samples*L,M) - pi;
    n=1:M;
    an = (theta + (2*pi*n - pi))/(4*M);
    
    %--Generating time trace--%
    
    
    t=0;
    
    %evaluate Xc(t)
    Xc = sqrt(2/M) * sum(cos((wd*t*cos(an)) + phi) .* cos(psi),2);

    %evaluate Xs(t)
    Xs = sqrt(2/M) * sum(cos((wd*t*cos(an)) + phi) .* sin(psi),2);

    X_0 = reshape(Xc + 1j*Xs,L,no_samples);
    
    t=k*Ts;
    
    %evaluate Xc(t)
    Xc = sqrt(2/M) * sum(cos((wd*t*cos(an)) + phi) .* cos(psi),2);

    %evaluate Xs(t)
    Xs = sqrt(2/M) * sum(cos((wd*t*cos(an)) + phi) .* sin(psi),2);

    X_k = reshape(Xc + 1j*Xs,L,no_samples);
    %size(X_k)
end