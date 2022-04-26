sdata = read_spec('eq.spec.h5')


%%

x1 = 0.3;
x2 = 0.4;
x3 = 0.5;
lvol = 1;

dx = 1e-8;

dx1 = x1 + dx; dx2 = x2 + dx; dx3 = x3 + dx;

"-------------NEW RUN-------------"

jac = get_spec_jacobian(sdata,lvol,x1,x2,x3);

djac = [get_spec_jacobian(sdata,lvol,dx1,x2,x3)-jac, ...
    get_spec_jacobian(sdata,lvol,x1,dx2,x3)-jac, ...
    get_spec_jacobian(sdata,lvol,x1,x2,dx3)-jac]/dx;

tmp = get_spec_magfield(sdata,lvol,x1,x2,x3);
for i = 1:3
    B(i) = tmp{i};
end

modB = get_spec_modB(sdata,lvol,x1,x2,x3);

tmp = get_spec_metric(sdata,lvol,x1,x2,x3);
for i = 1:3
    for j = 1:3
        met(i,j) = tmp{i}{j};
    end
end

tmp = get_spec_metric(sdata,lvol,dx1,x2,x3);
for i = 1:3
    for j = 1:3
        dGds(i,j) = (tmp{i}{j} - met(i,j))/dx;
    end
end
tmp = get_spec_metric(sdata,lvol,x1,dx2,x3);
for i = 1:3
    for j = 1:3
        dGdt(i,j) = (tmp{i}{j} - met(i,j))/dx;
    end
end
tmp = get_spec_metric(sdata,lvol,x1,x2,dx3);
for i = 1:3
    for j = 1:3
        dGtz(i,j) = (tmp{i}{j} - met(i,j))/dx;
    end
end




tmp = get_spec_magfield(sdata,lvol,dx1,x2,x3);
for i = 1:3
    dBds(i) = (tmp{i} - B(i))/dx;
end
tmp = get_spec_magfield(sdata,lvol,x1,dx2,x3);
for i = 1:3
    dBdt(i) = (tmp{i} - B(i))/dx;
end
tmp = get_spec_magfield(sdata,lvol,x1,x2,dx3);
for i = 1:3
    dBdz(i) = (tmp{i} - B(i))/dx;
end



jac
met
"metric derivatives"
dGds
dGdt
dGtz

B
modB

dBds
dBdt
dBdz



