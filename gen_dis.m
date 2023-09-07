function [d] = gen_dis(loc_a, loc_b)
d = 0;
for i=1:2 % ectangular coordinate system number of x, y
    d = d+(loc_a(i)-loc_b(i))^2;
end
d = sqrt(d);