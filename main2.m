function [] = main2(sch,L,modln,M,v,fc)
    close all;
    clc;
    tic
    L = str2double(L);
    M = str2double(M);
    fc = str2double(fc) * 1e6;
    v = str2double(v) *5/18;
    c=300000000;                    %speed of light     
    fd = fc*v/c;                    %Doppler Frequency
    Ts = 1/10*(1/(4*fd));
    
     
    k = 0 : 4 : 100 ;
   
    sample_size = 1e5;
    snr_db  = 15;
    snr = 10.^(snr_db/10);
    
    if modln == "MPSK"
        constellation = MPSK(M);
    elseif modln == "MQAM"
        constellation = MQAM(M);
        %sum(abs(constellation).^2)/length(constellation)
    end
    
    %unit energy samples
    
    
    Pe = zeros(1,length(snr_db));
    
    if sch == "MRC"
        x = datasample(constellation,sample_size);
        for i = 1 : length(k)
            [h_r,h_t] = Jake(v,fc,k(i),sample_size,L);
            n = sqrt(0.5)*randn(L,sample_size) + 1i*sqrt(0.5)*randn(L,sample_size); %CN(0,1) AWGN
            n = n/sqrt(snr);
            %L receiver system
            %PSK
            
            y = h_t .* x + n;
            %size(y)
            
            %finding column wise norm of h for normalization
            mod_hr = abs(h_r);
            norm2_hr = sum(mod_hr.^2,1);
            
            %MRC
            y_mrcx = conj(h_r) .* y;
            y_mrc = sum(y_mrcx,1);
            %find(isnan((y_mrc)))
            y_mrc = y_mrc./norm2_hr;
        
            if modln == "MPSK"
                x_hat = psk_dec(y_mrc, constellation, M);
            elseif modln == "MQAM"
                
                x_hat = MLD(y_mrc,constellation);
                
            end
            
            Pe(i) = sum(x_hat ~= x)/(sample_size);
            
        end
    
    
    elseif sch == "ALM"
        x1 = datasample(constellation,sample_size);
        x2 = datasample(constellation,sample_size);
        for i = 1 : length(k)
            [h1_r,h1_t] = Jake(v,fc,k(i),sample_size,L);
            [h2_r,h2_t] = Jake(v,fc,k(i),sample_size,L);
            
            n1 = sqrt(0.5)*randn(L,sample_size) + 1i*sqrt(0.5)*randn(L,sample_size); %CN(0,1) AWGN
            n1 = n1/sqrt(snr);

            n2 = sqrt(0.5)*randn(L,sample_size) + 1i*sqrt(0.5)*randn(L,sample_size); %CN(0,1) AWGN
            n2 = n2/sqrt(snr);

            mod_h1_r = abs(h1_r);
            norm2_h1_r = sum(mod_h1_r.^2,1);
            mod_h2_r = abs(h2_r);
            norm2_h2_r = sum(mod_h2_r.^2,1);
            
            %y1 time instant 1
            %y2 time instant 2
            y1 = h1_t.*x1 +h2_t .*x2 +n1;
            y2 = -h2_t.*conj(x1) + h1_t.*conj(x2) + n2;

            %decoding x1 and x2
            x1_ds = sum((conj(h1_r).*y1 - h2_r.*conj(y2)),1);
            x1_ds = x1_ds./(norm2_h2_r + norm2_h1_r);
            x2_ds = sum((conj(h2_r).*y1 + h1_r.*conj(y2)),1);
            x2_ds = x2_ds./(norm2_h2_r + norm2_h1_r);


            if modln == "MPSK"
                x1_hat = psk_dec(x1_ds, constellation, M);
                x2_hat = psk_dec(x2_ds, constellation, M);
            elseif modln == "MQAM"
                x1_hat = MLD(x1_ds,constellation);
                x2_hat = MLD(x2_ds,constellation);
                
            end
            Pe(i) = sum([(x1_hat ~=  x1),(x2_hat ~=  x2)])/(2*sample_size);
            
        end
    end
    figure
    semilogy(k*Ts*fd,Pe,'-r*')
    xlabel("kT_{s}f_{d}")
    ylabel("SEP")
    title([num2str(L),' Rx Antennas ',num2str(sch),' ',num2str(modln)])
    toc
    %QAM
end