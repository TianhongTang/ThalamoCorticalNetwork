function col = getUniqueColor(e,a)

main  = floor(e/4);
other = mod(e,4);
m = 0.5+main*0.12;
o1 = 0.25-other*0.06;
o2 = 0.25+other*0.06;

if a == 1
  col = [m o1 o2];
elseif a == 2
  col = [o1 m o2];
elseif a == 4
  col = [o1 o2 m];
else 
  col = [0 0 0];
end

