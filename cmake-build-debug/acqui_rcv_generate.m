nrcv    = 80;
drcv    = 24;
z0      = .5;
xrcv    = 0:drcv:drcv*nrcv; 
zrcv    = ones(length(xrcv),1) * z0; 

f = fopen("acqui_rcv", 'w');
for i = 1:length(xrcv)
    fprintf(f, '%.2f %.2f\n', zrcv(i),xrcv(i));
end

fclose(f);
