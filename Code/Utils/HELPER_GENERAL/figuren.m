function [h] = figuren(varargin)
% [h] = figuren(varargin)
% same as 'figure', but applies my favorite settings using 'nsubplot'!

h = figure(varargin{:});
nsubplot;

% custom settings added 2014-02-03

% put figure in upper-left corner
screensize = get(0,'screensize');
pos = get(h,'outerposition');
pos(1) = 1;
pos(2) = screensize(4)-pos(4);
set(h,'outerposition',pos);

% if ismac
%     % change the menu so that command-shift-s
%     %  activates "save as", like it does
%     %  in all other applications. We need to
%     %  do this because MATLAB for macs is
%     %  bizarre and only includes single-button
%     %  shortcuts in menus (e.g. command-s,
%     %  not command-shift-s). And even if you
%     %  try to edit the Accelerator keys for
%     %  the existing menus, it only lets you add
%     %  single-button accelerators.
%     %
%     % so we have to use this trick adapted from:
%     %  http://undocumentedmatlab.com/blog/customizing-menu-items-part-2/
%     %  to access the underlying Java objects
%     %  and set the accelerator key there.
%     jFrame = get(handle(h),'JavaFrame');
%     try
%         % R2008a and later
%         jMenuBar = jFrame.fHG1Client.getMenuBar;
%     catch
%         % R2007b and earlier
%         jMenuBar = jFrame.fFigureClient.getMenuBar;
%     end
%     
%     % File main menu is the first main menu item => index=0
%     jFileMenu = jMenuBar.getComponent(0);
% 
%     % Save menu item is the 6th menu item (separators included)
%     jSaveAs = jFileMenu.getMenuComponent(5); %Java indexes start with 0!
%     % just to be sure that this is really the 'Save As...' menu item
%     if isequal(get(jSaveAs,'Label'),'Save As...')
%         % set a new accelerator key for this menu item:
%         jAccelerator = javax.swing.KeyStroke.getKeyStroke('meta shift S');
%         jSaveAs.setAccelerator(jAccelerator);
%     else
%         warning('figuren on mac could not set the menu item "Save As..." to have the accelerator "meta shift S"');
%     end;
%     
%     % scale to be bigger, since I now have bigger screens.
%     scale(h,'scale',2);
% end;
