function constellation = MQAM(M)
    if log2(M) ~= floor(log2(M))
        constellation = -1;
        disp("M is not a power of 2")
        
    else
        r = sqrt(M);

        if r ~= floor(r)
            r = sqrt(2*M);
        end
       
        c = M/r;
        x = [0 : r-1]*2 - (r-1);
        y = [0 : c-1]*2 - (c-1);

        constellation =   x + 1i*y';    
        
        
        constellation = reshape(constellation.',1,[]);
        constellation = constellation * (sqrt(M)/sqrt(sum(abs(constellation).^2)));
        
        
    end
 
end