function dgz = loadDGZ(filename)
  
  [path,fileroot,ext] = fileparts(filename);
  DGZPath = 'f:/Data/Objects/dgz';
  SPKPath = 'f:/Data/Objects/spk';
  
  if isempty(path), path = DGZPath; end
  dgzfullpath = sprintf('%s/%s.dgz',path,fileroot);
  if isempty(dir(dgzfullpath))
    errordlg(sprintf('Cannot find file %s',dgzfullpath),'Load DGZ');
    return
  end
  dgz         = dg_read(dgzfullpath);

  spkfullpath = sprintf('%s/%s.spk',SPKPath,fileroot);  
  if ~isempty(dir(spkfullpath))
    dgz.spk     = spk_read(spkfullpath);
  end
  


%------------------
%
function spk = spk_read(file)
  spk = checkFile(file,'SPKT');
  