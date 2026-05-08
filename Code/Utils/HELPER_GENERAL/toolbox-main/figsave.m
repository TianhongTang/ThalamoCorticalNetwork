function figsave(filename)

% If we just want a quick save, make a random 16 alphanumeric string
if nargin < 1
    symbols = ['a':'z' 'A':'Z' '0':'9'];
    MAX_ST_LENGTH = 16;
    nums = randi(numel(symbols),[1 MAX_ST_LENGTH]);
    st = symbols (nums);
    saveas(gcf,['/Users/stevenerrington/Desktop/' st '.png'])

% Otherwise... if we have a filename that we want to save by, then use that 
else
    saveas(gcf,['/Users/stevenerrington/Desktop/' filename '.png'])
end

end
