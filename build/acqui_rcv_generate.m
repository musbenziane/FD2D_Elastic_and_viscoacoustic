nrcv    = 701;
drcv    = 4;
z0      = .5;
xrcv    = drcv*50:drcv:drcv*799-drcv*50;
zrcv    = ones(length(xrcv),1) * z0;

f = fopen("acqui_rcv", 'w');
for i = 1:length(xrcv)
fprintf(f, '%.2f %.2f\n', zrcv(i),xrcv(i));
end

fclose(f);