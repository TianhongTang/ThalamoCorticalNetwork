function acqp = PVrdPar(dir, filenum, optin);
% %  acqp = PVrdPar(dir, filenum, optin)
% %
% %  Bruker ParaVision
% %  reads parameter files ACQP, IMND from <dir,file>
% %
% %  routine is initialized by PVinitAcqp(): add any additional parameter to be read there !
% %
% %  default options:
% %     opt(  'VERBOSE','1')
% %
% %  return: struct acqp/imnd pars (flexible format, easy to extend)
% %
% %  Caveats: 
% %     2D par readout not implemented, e.g. ACQ_trim=(19, 3)
% %
% %  Sep 2001 -  Josef Pfeuffer
% %
FCTNAME = 'PVrdPar';

global STDPATH

PARCHAR = '##$';
PARFILE1 = 'acqp';     % needs path separator "\" upfront
PARFILE2 = 'imnd';
PARFILE3 = 'method';

      %%%%  option handling   %%%%%
      % --- default options: define whole struct
dopt.VERBOSE = 1;

      % --- arg handling
nargVars = 2;
narg = nargin;
error(nargchk(nargVars,nargVars+1,narg))
if (narg == nargVars+1)
   dopt = setopt(dopt,optin);
   narg = narg - 1;    % narg is NOT including options
end
if (filenum > 0)    % stdpath handling
   file1 = sprintf('%s%s/%s/%s', STDPATH.pv,num2str(dir),num2str(filenum),PARFILE1);
   file2 = sprintf('%s%s/%s/%s', STDPATH.pv,num2str(dir),num2str(filenum),PARFILE2);
   file3 = sprintf('%s%s/%s/%s', STDPATH.pv,num2str(dir),num2str(filenum),PARFILE3);
else              % give full path as arg1
   file1 = sprintf('%s/%s', num2str(filenum), PARFILE1);
   file2 = sprintf('%s/%s', num2str(filenum), PARFILE2);
   file3 = sprintf('%s/%s', num2str(filenum), PARFILE3);
end
fname = {file1,file2,file3};
numParFiles = length(fname);

f_verbose = dopt.VERBOSE;
      %%%%  end: handling options  %%%%%

acqp      = PVinitAcqp;
acqpNames = fieldnames(acqp);
acqpNtags = length(acqpNames);

for ifile=1:numParFiles
   
filename = char(fname(ifile));
fid = fopen(filename,'r');
if (fid > 0)
   if f_verbose fprintf('reading <%s>\n', filename); end
   
  while (~feof(fid))
   tline = fgetl(fid);
   %%disp(tline)
   strncmp(PARCHAR, tline, length(PARCHAR));
   if (strncmp(PARCHAR, tline, length(PARCHAR)))   % check if parameter
      par = tline( length(PARCHAR)+1 : findstr(tline,'=')-1 );
      parvalStr = tline( findstr(tline,'=')+1 : length(tline) );
      for itag=1:acqpNtags      % searching acqp STRUCT
         if (strcmp(par, acqpNames(itag)))      
                           % discriminate string/ numeric
                           % --- string handling 
            if ( isstr( getfield(acqp, char(acqpNames(itag))) ) )
               if isempty( findstr('(', parvalStr) )  
                  parval = parvalStr;
               else        % read only until the first '>'
                  parval = '';                 
                  ch = fscanf(fid, '%c',1);
                  while (ch ~= '#')
                     parval = [parval ch];
                     if (ch == '>')
                        tline = fgetl(fid);     % get rest of line
                     end
                     ch = fscanf(fid, '%c',1);
                  end
                  stat = fseek(fid,-1,'cof');   % reposition fid
               end              
            else           % --- numeric handling 
               if isempty( findstr('(', parvalStr) )  
                  parval = str2num(parvalStr);
               else        % array handling, read formatted value(s)
                  parval = fscanf(fid, '%f', numeric(parvalStr));
               end
            end
            acqp = setfield(acqp, char(acqpNames(itag)), parval);
            %%disp(par)
            %%disp(parvalStr)
            %%whos parval
         end   %
      end   % for 
   end   % else ignore line
  end
  fclose(fid);

else
   if f_verbose fprintf('can not open <%s>\n', filename); end
end


end %% for ifile
%%parval = cellstr(parval);

%------------------------------------------------------------
function acqp = PVinitAcqp;
% %  acqp = PVinitAcqp
% %
% %  Bruker ParaVision
% %  inits acqp struct and defines thereby the ACQP,IMND,METHOD parameter to be read
% %
% %  add any parameter from ACQP,IMND,METHOD
% %
% %
FCTNAME = 'PVinitAcqp';

   % define struct to be read here:
   %     = 0   type: double
   %     = ''  type: string
   %           ( arrays are determined automatically from '=( x )' => next lines interpreted )
acqp.PULPROG = '';
acqp.GRDPROG = '';
acqp.ACQ_dim = 0;
acqp.ACQ_dim_desc = '';
acqp.ACQ_size = 0;
acqp.ACQ_word_size = '';
acqp.BYTORDA = '';
acqp.NSLICES = 0;
acqp.ACQ_ns_list_size = 0;
acqp.ACQ_ns = 0;
acqp.ACQ_ns_list = 0;
acqp.NECHOES = 0;
acqp.NS = 0;
acqp.NI = 0;
acqp.NA = 0;
acqp.NAE = 0;
acqp.NR = 0;
acqp.DS = 0;
acqp.D = 0;
acqp.L = 0;
acqp.SP = 0;
acqp.RG = 0;
acqp.SW_h = 0;
acqp.SFO1 = 0;
acqp.DIGMOD = '';
acqp.DECIM = 0;
acqp.DSPFVS = 0;
acqp.ACQ_vd_list = 0;
acqp.IMND_matrix = 0;
acqp.IMND_n_slices = 0;
acqp.IMND_n_echo_images = 0;
acqp.IMND_echo_scan_mode = '';
acqp.IMND_echo_scan_eq = 0;
acqp.IMND_echo_time = 0;
acqp.IMND_n_echos = 0;
acqp.IMND_n_averages = 0;
acqp.IMND_rep_time = 0;
acqp.IMND_recov_time = 0;
acqp.IMND_pulse_length = 0;
acqp.EPI_resid = '';
acqp.EPI_segmentation_mode = '';
acqp.IMND_numsegments = 0;
acqp.EPI_scan_mode = '';
acqp.EPI_TE_eff = 0;
acqp.EPI_zero_phase_ms = 0;
acqp.EPI_seg_acq_time = 0;
acqp.EPI_slice_rep_time = 0;
acqp.EPI_image_rep_time = 0;
acqp.EPI_image_time = 0;
acqp.EPI_swh_eff_phase = 0;
acqp.EPI_use_vd = '';
acqp.EPI_use_id = '';
acqp.EPI_vd_start = 0;
acqp.IMND_slice_thick = 0;
acqp.IMND_slicepack_position = 0;
acqp.IMND_slicepack_gap = 0;
acqp.IMND_fov = 0;
acqp.IMND_sl_thick_hz = 0;
acqp.IMND_acq_time = 0;
acqp.IMND_DW_time = 0;
acqp.MP_Mode = '';
acqp.PVM_EchoTime = 0;
acqp.MP_MixingTime = 0;
acqp.PVM_RepetitionTime = 0;
acqp.PVM_NAverages = 0;
acqp.PVM_NRepetitions = 0;
acqp.PVM_NEchoImages = 0;
acqp.EffectiveTE = 0;
acqp.MP_AcqSize = 0;
acqp.PVM_EffSWh = 0;
acqp.PVM_SliceThick = 0;
% --- end of acqp/imnd field definition !
