function z = util_logsumexp(x)
% z = log(sum(exp(x)))

x=sort(x,'descend');
if x(1) == -Inf
   z = -Inf;
else
   z = x(1)+log(1+sum(exp(x(2:end)-x(1))));
end

