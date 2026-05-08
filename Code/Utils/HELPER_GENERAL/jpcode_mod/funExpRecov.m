function F = funExpRecov(P,x)

if length(P) == 2
    F  = P(1).*(1 - exp(-x./P(2)));
elseif length(P) == 3
    F  = P(1).*(1 - P(3)*exp(-x./P(2)));
else
    stop
end