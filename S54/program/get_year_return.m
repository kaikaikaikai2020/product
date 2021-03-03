
function temp = get_year_return(tref,ret) 
y=year(datenum(tref));
y_u = unique(y);
temp = zeros(size(y_u));
for i = 1:length(y_u)
    temp1 = ret(eq(y,y_u(i)));
    temp1(1) = 0;
    temp1 = cumprod(1+temp1);
    temp(i) = (temp1(end)-1)*100;
end
temp = num2cell([y_u,y_u,temp]);
end