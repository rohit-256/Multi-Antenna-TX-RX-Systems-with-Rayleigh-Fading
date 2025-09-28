%%%%%%% one transmit L receive %%%%%%

function [] = main1(sch,L,modln,M)
clc;
close all;
   
    tic

   
    sample_size = 1e6;
    snr_db  = -20:1:50;
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
        for i = 1 : length(snr)
            h = sqrt(0.5)*randn(L,sample_size) + 1i*sqrt(0.5)*randn(L,sample_size); %CN(0,1)
            n = sqrt(0.5)*randn(L,sample_size) + 1i*sqrt(0.5)*randn(L,sample_size); %CN(0,1) AWGN
            n = n/sqrt(snr(i));
            %L receiver system
            %PSK
            
            y = h .* x + n;
            %size(y)
            
            %finding column wise norm of h for normalization
            mod_h = abs(h);
            norm2_h = sum(mod_h.^2,1);
            
            %MRC
            y_mrcx = conj(h) .* y;
            y_mrc = sum(y_mrcx,1);
            %find(isnan((y_mrc)))
            
        
            if modln == "MPSK"
                x_hat = psk_dec(y_mrc, constellation, M);
            elseif modln == "MQAM"
                y_mrc = y_mrc./norm2_h;
                x_hat = MLD(y_mrc,constellation);
                
            end
            
            Pe(i) = sum(x_hat ~= x)/(sample_size);
            if(Pe(i) == 0)
                break;
            end
        end
    
    
    elseif sch == "ALM"
        x1 = datasample(constellation,sample_size);
        x2 = datasample(constellation,sample_size);
        for i = 1 : length(snr)
            h1 = sqrt(0.5)*randn(L,sample_size) + 1i*sqrt(0.5)*randn(L,sample_size); %CN(0,1)
            mod_h1 = abs(h1);
            norm2_h1 = sum(mod_h1.^2,1);

            h2 = sqrt(0.5)*randn(L,sample_size) + 1i*sqrt(0.5)*randn(L,sample_size); %CN(0,1)
            mod_h2 = abs(h2);
            norm2_h2 = sum(mod_h2.^2,1);
            n1 = sqrt(0.5)*randn(L,sample_size) + 1i*sqrt(0.5)*randn(L,sample_size); %CN(0,1) AWGN
            n1 = n1/sqrt(snr(i));

            n2 = sqrt(0.5)*randn(L,sample_size) + 1i*sqrt(0.5)*randn(L,sample_size); %CN(0,1) AWGN
            n2 = n2/sqrt(snr(i));

            
            
            %y1 time instant 1
            %y2 time instant 2
            y1 = h1.*x1 +h2 .*x2 +n1;
            y2 = -h2.*conj(x1) + h1.*conj(x2) + n2;

            %decoding x1 and x2
            x1_ds = sum((conj(h1).*y1 - h2.*conj(y2)),1);
            x1_ds = x1_ds./(norm2_h2 + norm2_h1);
            x2_ds = sum((conj(h2).*y1 + h1.*conj(y2)),1);
            x2_ds = x2_ds./(norm2_h2 + norm2_h1);


            if modln == "MPSK"
                x1_hat = psk_dec(x1_ds, constellation, M);
                x2_hat = psk_dec(x2_ds, constellation, M);
            elseif modln == "MQAM"
                x1_hat = MLD(x1_ds,constellation);
                x2_hat = MLD(x2_ds,constellation);
                
            end
            Pe(i) = sum([(x1_hat ~=  x1),(x2_hat ~=  x2)])/(2*sample_size);
            if(Pe(i) == 0)
                break;
            end
        end
    end
    figure
    semilogy(snr_db,Pe,'-r*')
    xlabel("SNR")
    ylabel("SEP")
    title([num2str(L),' Rx Antennas ',num2str(sch),' ',num2str(modln)])
    toc
end
