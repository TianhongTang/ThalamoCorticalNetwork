function events()
%
% This function sets EVT as a global structure
%

global EVT

EVT.BEGINOBS        =             19;
EVT.ENDOBS          =             20;
EVT.TRIALTYPE       =             22;
EVT.OBSTYPE         =             23;
EVT.EMLOG           =             24;
EVT.FIXSPOT         =             25;
EVT.EMPARAMS        =             26;
EVT.STIMULUS        =             27;
EVT.PATTERN         =             28;
EVT.STIMTYPE        =             29;
EVT.CUE		    =             32;
EVT.SOUND	    =             35;
EVT.FIXATE	    =             36;
EVT.RESPONSE        =             37;
EVT.ABORT	    = 		  41;
EVT.REWARD	    = 		  42;
EVT.PUNISH	    = 		  44;
EVT.PHYS            =   45;
% added 11-0703 AM
EVT.MRI		    = 		  46;
EVT.FLOATS_1        =            128;
EVT.FLOATS_2        =            129;
EVT.FLOATS_3        =            130;
EVT.FLOATS_4        =            131;
EVT.FLOATS_5        =            132;
EVT.FLOATS_6        =            200;
EVT.FLOATS_7        =            201;
EVT.FLOATS_8        =            202;
EVT.FLOATS_9        =            203;
EVT.FLOATS_10        =           204;
EVT.FLOATS_11        =           205;
EVT.FLOATS_12        =           206;
EVT.FLOATS_13        =           207;
EVT.FLOATS_14        =           208;
EVT.FLOATS_15        =           209;
EVT.STRINGS_1       =            133;
EVT.STRINGS_2       =            134;
EVT.STRINGS_3       =            135;
EVT.RESETDISPLAY    =            133;
EVT.PATTERN_ON      =              1;
EVT.PATTERN_OFF     =              0;

