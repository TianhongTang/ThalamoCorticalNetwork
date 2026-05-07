function si2 = sortSession(si)

[a,indx] = sort(si.areas);
si2.electrodes = si.electrodes(indx);
si2.areas      = a;
