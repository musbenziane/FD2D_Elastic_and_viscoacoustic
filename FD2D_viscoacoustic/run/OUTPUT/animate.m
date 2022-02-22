clear; clc; close all

f=fopen("field_P_00001","r");
u=fread(f,"float64");
fclose(f);

u = reshape(u,260,[],320);

for i=1:120
 imagesc(reshape(u(i,:,:),[],320));
 pause(.09); 
 end