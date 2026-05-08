function srtlist = sortByLists(data, list1, list2, list3, list4)

if (nargin < 2)
  error('usage: srtlist = sortByLists(data, list1, [list2, list3, list4]');
end

sorted = [];
if nargin == 2
   tmp = unique(list1);
   for i = 1:length(tmp)
	sorted{i} = data(list1 == tmp(i));
   end	  	   
elseif nargin == 3
   tmp1 = unique(list1);
   tmp2 = unique(list2);
   for j = 1:length(tmp2)
      for i = 1:length(tmp1)
	sorted{(j-1)*length(tmp1)+i}...
		 = data(list1 == tmp1(i) & list2 == tmp2(j));
      end
   end	  	   
elseif nargin == 4
   tmp1 = unique(list1);
   tmp2 = unique(list2);
   tmp3 = unique(list3);
   for k = 1:length(tmp3)
      for j = 1:length(tmp2)
         for i = 1:length(tmp1)
	   sorted{(k-1)*length(tmp2)*length(tmp1)+(j-1)*length(tmp1)+i}...
  	      = data(list1 == tmp1(i) & list2 == tmp2(j) & list3 == tmp3(k));
         end
      end	  	   
   end	
end

srtlist = sorted;


