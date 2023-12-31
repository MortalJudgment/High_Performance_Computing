%DEMO2  MGLab Fourier Analysis Demo
%
%       Poisson problem on a 49x49 mesh.
%       The initial guess has a mix of low and high
%       frequencies.
%
%       V-Cycle parameters:
%
%       Two levels (two grid algorithm)
%       Gauss-Seidel smoothing; nu1 = 0;  nu2 = 4;
%       Half-weighting; Cubic interpolation
%       Sparse Gaussian elimination on coarse grid
%
%       We examine the error in Physical and Fourier
%       space two times in each V-Cycle: After the
%       coarse grid correction and after the post-
%       smoothings.
%
%       Accesses global variables in "demo_globals"

% James Bordner and Faisal Saied
% Department of Computer Science
% University of Illinois at Urbana-Champaign
% 10 April 1995

function demo2_uic_continue = demo2_info

demo_globals

page=[	'                                                 ';...
	'                    Demo 2                       ';...
        '                                                 ';...
	' Poisson problem on a 49x49 mesh.                ';...
	' The initial guess has a mix of low and high     ';...
	' frequencies.                                    ';...
	'                                                 ';...
	' V-Cycle parameters:                             ';...
	'                                                 ';...
	' Two levels (two grid algorithm)                 ';...
	' Gauss-Seidel smoothing; nu1 = 0;  nu2 = 4;      ';...
	' Half-weighting; Cubic interpolation             ';...
	' Sparse Gaussian elimination on coarse grid      ';...
	'                                                 ';...
	' We examine the error in Physical and Fourier    ';...
	' space two times in each V-Cycle: After the      ';...
	' coarse grid correction and after the post-      ';...
	' smoothings.                                     ';...
	'                                                 ';...
        '                                                 ';...
];

DEMO_INFO_FIG = figure('Position', [30 30 500 600],...
	   'Name', 'Demo2 Info',...
	   'NumberTitle', 'off', ...
	   'Color','blue');


[mm,nn] = size(page);

subplot(1,1,1)
hold off
cla

ht = 0.95;
for j = 1:mm
   text(0.05, ht, page(j,:))
   axis('off')
   ht = ht - 0.04;
end

cb_uic_continue = 'close(DEMO_INFO_FIG), demo2_run';

demo2_uic_continue = uicontrol(...
		'units', 'normalized',...
        	'Position', [0.1 0.15 0.35 0.07],...
		'style', 'pushbutton', ...
		'string', 'Click here to continue',...
		'backgroundcolor', 'cyan',...
		'foregroundcolor', [0 0 0],...
		'callback', cb_uic_continue);

cb_uic_cancel = 'close(DEMO_INFO_FIG)';

demo1_uic_cancel = uicontrol(	...
		'units', 'normalized',...
        'Position', [0.55 0.15 0.35 0.07],...
		'style', 'pushbutton', ...
		'string', 'Cancel',...
		'backgroundcolor', 'cyan',...
		'foregroundcolor', [0 0 0],...
		'callback', cb_uic_cancel);



