%MGLab V0.00beta   Interactive Multigrid Package

% James Bordner and Faisal Saied
% Department of Computer Science
% University of Illinois at Urbana-Champaign
% 10 April 1995

include_flags
include_globals
include_figs
demo_globals

% Initialize parameter defaults

set_defaults;

% == MAIN MENU =================================================

bgc = [0.9 0.9 1.0];

main_fig = figure('Position', main_position,...
    'Name', 'MGLab',...
    'NumberTitle', 'off', ...
    'Color','yellow');

% == MGLab Menu Item ===========================================

f_mglab = menu_header(main_fig,'MGLab','on','on');
    menu_item(f_mglab,'Run', 'off','on','[sol1,resids1,its1]=run;');
    menu_item(f_mglab,'Show Params','off','on','show_params;');
    menu_item(f_mglab,'Version Info','off','on','version_info;');
    menu_item(f_mglab,'Reset','off','on','set_defaults;');
    menu_item(f_mglab,'Restart','off','on','close(main_fig); close; MGLab');
    menu_item(f_mglab,'Quit','off','on','close(main_fig); close');

% == Problem Menu Item ============================================

f_problem=menu_header(main_fig,'Problem','on','on');

    menu_item(f_problem,'Poisson','on','on',...
        'problem_flag = POISSON;generate_matrix=1;');
    f_problem_1 = menu_item(f_problem,'Helmholtz', 'off','on',...
        'problem_flag = HELMHOLTZ;generate_matrix=1;prob_args(1) = -10;');
    f_problem_11=menu_header(f_problem_1,'k = ','on','on');
        menu_item(f_problem_11,'-10','off','on',...
            'prob_args(1)=-10;');
        menu_item(f_problem_11,'-5','off','on',...
            'prob_args(1)=-5;');
        menu_item(f_problem_11,'-1','off','on',...
            'prob_args(1)=-1;');
        menu_item(f_problem_11,'0','off','on',...
            'prob_args(1)=0;');
        menu_item(f_problem_11,'1','off','on',...
            'prob_args(1)=1;');
        menu_item(f_problem_11,'5','off','on',...
            'prob_args(1)=5;');
        menu_item(f_problem_11,'10','off','on',...
            'prob_args(1)=10;');
        menu_item(f_problem_11,'10+ i','off','on',...
            'prob_args(1)=10+sqrt(-1);');
    f_problem_2 = menu_item(f_problem,'Convection-Diffusion', 'off','on',...
    'problem_flag=CONVECT_DIFFUSE;generate_matrix=1;');
    f_problem_21=menu_header(f_problem_2,'Lambda = ','on','on');
        menu_item(f_problem_21,'0','off','on',...
            'prob_args(1)=0;');
        menu_item(f_problem_21,'10','off','on',...
            'prob_args(1)=10;');
        menu_item(f_problem_21,'100','off','on',...
            'prob_args(1)=100;');
        menu_item(f_problem_21,'1000','off','on',...
            'prob_args(1)=1000;');
    f_problem_22=menu_header(f_problem_2,'Sigma = ','on','on');
        menu_item(f_problem_22,'0','on','on',...
            'prob_args(2)=0;');
        menu_item(f_problem_22,'5','off','on',...
            'prob_args(2)=5;');
        menu_item(f_problem_22,'10','off','on',...
            'prob_args(2)=10;');
        menu_item(f_problem_22,'20','off','on',...
            'prob_args(2)=20;');
        menu_item(f_problem_22,'50','off','on',...
            'prob_args(2)=50;');
        menu_item(f_problem_22,'100','off','on',...
            'prob_args(2)=100;');
        menu_item(f_problem_22,'-50','off','on',...
            'prob_args(2)=-50;');
        menu_item(f_problem_22,'-100','off','on',...
            'prob_args(2)=-100;');
    f_problem_3=menu_item(f_problem,'Cut Square', 'off','on',...
    'problem_flag = CUT_SQUARE;generate_matrix=1;prob_args(1) = 10;');
    f_problem_31=menu_header(f_problem_3,'Alpha = ','on','on');
        menu_item(f_problem_31,'0.001','off','on',...
            'prob_args(1)=0.001;');
        menu_item(f_problem_31,'0.01','off','on',...
            'prob_args(1)=0.01;');
        menu_item(f_problem_31,'0.1','off','on',...
            'prob_args(1)=0.1;');
        menu_item(f_problem_31,'1','off','on',...
            'prob_args(1)=1;');
        menu_item(f_problem_31,'10','off','on',...
            'prob_args(1)=10;');
        menu_item(f_problem_31,'100','off','on',...
            'prob_args(1)=100;');
        menu_item(f_problem_31,'1000','off','on',...
            'prob_args(1)=1000;');
    menu_item(f_problem,'Poisson-Boltzmann', 'off','off',...
        'problem_flag=POISSON_BOLTZMAN;generate_matrix=1;');
    f_problem_4=menu_header(f_problem,'Problem Size','off','on');
        menu_item(f_problem_4,'  7   ','off','on',...
            [['nx1=7;ny1=7;generate_matrix=1;generate_rhs=1;']';...
            ['coarse_level=min([coarse_level max_level(nx1)]);']']');
        menu_item(f_problem_4,' 15   ','off','on',...
            [['nx1=15;ny1=15;generate_matrix=1;generate_rhs=1;']';...
            ['coarse_level=min([coarse_level max_level(nx1)]);']']');
        menu_item(f_problem_4,' 31   ','off','on',...
            [['nx1=31;ny1=31;generate_matrix=1;generate_rhs=1;']';...
            ['coarse_level=min([coarse_level max_level(nx1)]);']']');
        menu_item(f_problem_4,' 63   ','off','on',...
            [['nx1=63;ny1=63;generate_matrix=1;generate_rhs=1;']';...
            ['coarse_level=min([coarse_level max_level(nx1)]);']']');
        menu_item(f_problem_4,'127   ','off','on',...
            [['nx1=127;ny1=127;generate_matrix=1;generate_rhs=1;']';...
            ['coarse_level=min([coarse_level max_level(nx1)]);']']');
        menu_item(f_problem_4,'255   ','off','on',...
            [['nx1=255;ny1=255;generate_matrix=1;generate_rhs=1;']';...
            ['coarse_level=min([coarse_level max_level(nx1)]);']']');

% == Solver Menu Item ==============================================

f_solver = menu_header(main_fig,'Solver','on','on');

    menu_item(f_solver,'V-Cycle','off','on',...
        'solver_flag = C_VMG;');
    menu_item(f_solver,'PCG','off','on',...
        'solver_flag = C_PCG;');
    menu_item(f_solver,'BiCG-STAB','off','on',...
        'solver_flag = C_BICG_STAB;');
    menu_item(f_solver,'CGS','off','on',...
        'solver_flag = CGS;');
    menu_item(f_solver,'TFQMR','off','off',...
        'solver_flag = TFQMR;');

f_solver_1=menu_item(f_solver,'GMRES(k)','off','on',...
    'solver_flag = C_GMRES;');
f_solver_11=menu_header(f_solver_1,'k = ','on','on');
    menu_item(f_solver_11,'1','off','on','restart=1;');
    menu_item(f_solver_11,'5','off','on','restart=5;');
    menu_item(f_solver_11,'10','off','on','restart=10;');
    menu_item(f_solver_11,'15','off','on','restart=15;');
    menu_item(f_solver_11,'20','off','on','restart=20;');

f_solver_2 = menu_item(f_solver,'SOR','off','on',...
    'solver_flag = C_SOR;');
f_solver_21=menu_header(f_solver_2,'omega = ','on','on');
    menu_item(f_solver_21,'1','off','on','SOR_omega=1;');
    menu_item(f_solver_21,'1.1','off','on','SOR_omega=1.1;');
    menu_item(f_solver_21,'1.2','off','on','SOR_omega=1.2;');
    menu_item(f_solver_21,'1.3','off','on','SOR_omega=1.3;');
    menu_item(f_solver_21,'1.4','off','on','SOR_omega=1.4;');
    menu_item(f_solver_21,'1.5','off','on','SOR_omega=1.5;');
    menu_item(f_solver_21,'1.6','off','on','SOR_omega=1.6;');
    menu_item(f_solver_21,'1.7','off','on','SOR_omega=1.7;');
    menu_item(f_solver_21,'1.8','off','on','SOR_omega=1.8;');
    menu_item(f_solver_21,'1.9','off','on','SOR_omega=1.9;');

menu_item(f_solver,'Full-Multigrid','on','on',...
    'solver_flag = C_FMG;');

f_solver_precon=menu_header(f_solver,'Preconditioner','on','on');
    menu_item(f_solver_precon,'V-Cycle','off','on',...
        'precon_flag = MG_CYCLE;');
    menu_item(f_solver_precon,'Jacobi','off','on',...
        'precon_flag = JACOBI;');
    menu_item(f_solver_precon,'Block-Jacobi','off','off',...
        'precon_flag = BLOCK_JACOBI;');
    menu_item(f_solver_precon,'Gauss-Seidel','off','on',...
        'precon_flag = GAUSS_SEIDEL;');
    menu_item(f_solver_precon,'ILU','off','off',...
        'precon_flag = ILU');
    menu_item(f_solver_precon,'SSOR','off','off',...
        'precon_flag = SSOR');
    menu_item(f_solver_precon,'None','off','on',...
        'precon_flag = NONE;');

f_solver_stop=menu_header(f_solver,'Stopping Criteria','off','on');
f_stop_1=menu_header(f_solver_stop,'Residual Tolerance','on','off');
    menu_item(f_stop_1,'None','off','on','rtol=0;');
    menu_item(f_stop_1,'1e-1','off','on','rtol=1e-1;');
    menu_item(f_stop_1,'1e-2','off','on','rtol=1e-2;');
    menu_item(f_stop_1,'1e-3','off','on','rtol=1e-3;');
    menu_item(f_stop_1,'1e-4','off','on','rtol=1e-4;');
    menu_item(f_stop_1,'1e-5','off','on','rtol=1e-5;');
    menu_item(f_stop_1,'1e-6','off','on','rtol=1e-6;');
    menu_item(f_stop_1,'1e-7','off','on','rtol=1e-7;');
    menu_item(f_stop_1,'1e-8','off','on','rtol=1e-8;');
    menu_item(f_stop_1,'1e-9','off','on','rtol=1e-9;');
    menu_item(f_stop_1,'1e-10','off','on','rtol=1e-10;');
    menu_item(f_stop_1,'1e-12','off','on','rtol=1e-12;');
    menu_item(f_stop_1,'1e-14','off','on','rtol=1e-14;');
    menu_item(f_stop_1,'1e-16','off','on','rtol=1e-16;');
f_stop_2=menu_header(f_solver_stop,'(Precon) Residual Tolerance',...
    'off','on');
    menu_item(f_stop_2,'None','off','on','prtol=0;');
    menu_item(f_stop_2,'1e-1','off','on','prtol=1e-1;');
    menu_item(f_stop_2,'1e-2','off','on','prtol=1e-2;');
    menu_item(f_stop_2,'1e-3','off','on','prtol=1e-3;');
    menu_item(f_stop_2,'1e-4','off','on','prtol=1e-4;');
    menu_item(f_stop_2,'1e-5','off','on','prtol=1e-5;');
    menu_item(f_stop_2,'1e-6','off','on','prtol=1e-6;');
    menu_item(f_stop_2,'1e-7','off','on','prtol=1e-7;');
    menu_item(f_stop_2,'1e-8','off','on','prtol=1e-8;');
    menu_item(f_stop_2,'1e-9','off','on','prtol=1e-9;');
    menu_item(f_stop_2,'1e-10','off','on','prtol=1e-10;');
    menu_item(f_stop_2,'1e-12','off','on','prtol=1e-12;');
    menu_item(f_stop_2,'1e-14','off','on','prtol=1e-14;');
    menu_item(f_stop_2,'1e-16','off','on','prtol=1e-16;');

f_stop_3=menu_header(f_solver_stop,'Iteration Limit','off','on');
    menu_item(f_stop_3,' None','off','on','max_it=0;');
    menu_item(f_stop_3,'    1','off','on','max_it=1;');
    menu_item(f_stop_3,'    2','off','on','max_it=2;');
    menu_item(f_stop_3,'    3','off','on','max_it=3;');
    menu_item(f_stop_3,'    5','off','on','max_it=5;');
    menu_item(f_stop_3,'   10','off','on','max_it=10;');
    menu_item(f_stop_3,'   20','off','on','max_it=20;');
    menu_item(f_stop_3,'   30','off','on','max_it=30;');
    menu_item(f_stop_3,'   50','off','on','max_it=50;');
    menu_item(f_stop_3,'  100','off','on','max_it=100;');
    menu_item(f_stop_3,'  200','off','on','max_it=200;');
    menu_item(f_stop_3,'  300','off','on','max_it=300;');
    menu_item(f_stop_3,'  500','off','on','max_it=500;');
    menu_item(f_stop_3,' 1000','off','on','max_it=1000;');
f_stop_4=menu_header(f_solver_stop,'Time Limit','off','off');
    menu_item(f_stop_4,'None','off','on','max_time=0;');
    menu_item(f_stop_4,'1 sec','off','on','max_time=1;');
    menu_item(f_stop_4,'5 sec','off','on','max_time=5;');
    menu_item(f_stop_4,'10 sec','off','on','max_time=10;');
    menu_item(f_stop_4,'30 sec','off','on','max_time=30;');
    menu_item(f_stop_4,'1 min','off','on','max_time=1*60;');
    menu_item(f_stop_4,'5 min','off','on','max_time=5*60;');
    menu_item(f_stop_4,'10 min','off','on','max_time=10*60;');
    menu_item(f_stop_4,'30 min','off','on','max_time=30*60;');
    menu_item(f_stop_4,'1 hour','off','on','max_time=60*60;');
f_stop_5=menu_header(f_solver_stop,'MFlop Limit','off','off');
    menu_item(f_stop_5,'None','off','on','max_mflop=0;');
    menu_item(f_stop_5,'   1','off','on','max_mflop=1;');
    menu_item(f_stop_5,'   5','off','on','max_mflop=5;');
    menu_item(f_stop_5,'  10','off','on','max_mflop=10;');
    menu_item(f_stop_5,'  20','off','on','max_mflop=20;');
    menu_item(f_stop_5,'  50','off','on','max_mflop=50;');
    menu_item(f_stop_5,' 100','off','on','max_mflop=100;');

% == MG Parameters ===========================================

f_solver_mg=menu_header(main_fig,'MG-Parameters','on','on');
f_mg_1 = menu_header(f_solver_mg,'Number of Levels','on','on');
    menu_item(f_mg_1,'1','off','on',...
        'coarse_level=min([1,max_level(nx1)]); generate_matrix=1;');
    menu_item(f_mg_1,'2','off','on',...
        'coarse_level=min([2,max_level(nx1)]); generate_matrix=1;');
    menu_item(f_mg_1,'3','off','on',...
        'coarse_level=min([3,max_level(nx1)]); generate_matrix=1;');
    menu_item(f_mg_1,'4','off','on',...
        'coarse_level=min([4,max_level(nx1)]); generate_matrix=1;');
    menu_item(f_mg_1,'5','off','on',...
        'coarse_level=min([5,max_level(nx1)]); generate_matrix=1;');
f_mg_2=menu_header(f_solver_mg,'Smoother','off','on');
f_mg_21=menu_item(f_mg_2,'Weighted Jacobi','on','on',...
    'smooth_flag=WEIGHTED_JACOBI;');
f_mg_211=menu_header(f_mg_21,'Weight = ','on','on');
    menu_item(f_mg_211,'1.00','off','on','wt=1.0;');
    menu_item(f_mg_211,'0.95','off','on','wt=0.95;');
    menu_item(f_mg_211,'0.90','off','on','wt=0.90;');
    menu_item(f_mg_211,'0.85','off','on','wt=0.85;');
    menu_item(f_mg_211,'0.80','off','on','wt=0.80;');
    menu_item(f_mg_2, 'Gauss-Seidel','off','on',...
    'smooth_flag=GAUSS_SEIDEL;');
    menu_item(f_mg_2, 'Red/Black Gauss-Seidel','off','off',...
        'smooth_flag=RB_GAUSS_SEIDEL;');
f_mg_22=menu_header(f_mg_2,'Pre-smoothings','on','on');
    menu_item(f_mg_22,'0','off','on','nu1=0;');
    menu_item(f_mg_22,'1','off','on','nu1=1;');
    menu_item(f_mg_22,'2','off','on','nu1=2;');
    menu_item(f_mg_22,'3','off','on','nu1=3;');
    menu_item(f_mg_22,'4','off','on','nu1=4;');
    menu_item(f_mg_22,'5','off','on','nu1=5;');
f_mg_23=menu_header(f_mg_2,'Post-smoothings','off','on');
    menu_item(f_mg_23,'0','off','on','nu2=0;');
    menu_item(f_mg_23,'1','off','on','nu2=1;');
    menu_item(f_mg_23,'2','off','on','nu2=2;');
    menu_item(f_mg_23,'3','off','on','nu2=3;');
    menu_item(f_mg_23,'4','off','on','nu2=4;');
    menu_item(f_mg_23,'5','off','on','nu2=5;');

f_mg_3=menu_header(f_solver_mg,'Restriction','off','on');
    menu_item(f_mg_3, 'Injection','off','on',...
        'restrict_flag=INJECTION;');
    menu_item(f_mg_3, 'Half Weighting','off','on',...
        'restrict_flag=HALF_WEIGHTING;');
    menu_item(f_mg_3, 'Full Weighting','off','on',...
        'restrict_flag=FULL_WEIGHTING;');
    menu_item(f_mg_3, 'Bilinear Adjoint','off','off',...
        'restrict_flag=BILINEAR_ADJOINT;');

f_mg_4=menu_header(f_solver_mg,'Prolongation','off','on');
    menu_item(f_mg_4, 'Linear','off','on',...
        'interp_flag=LINEAR;');
    menu_item(f_mg_4, 'Cubic','off','on',...
        'interp_flag=C_CUBIC;');
    menu_item(f_mg_4, 'Operator-based','off','off',...
        'interp_flag=OPERATOR_BASED;');
    menu_item(f_mg_4, 'Explicit/Bilinear','off','off',...
        'interp_flag=EXPLICIT_BILINEAR;');

f_mg_5=menu_header(f_solver_mg,'Coarse-grid Solver','off','on');
    menu_item(f_mg_5,'Sparse GE','off','on',...
        'coarse_solver_flag=DIRECT;');
    menu_item(f_mg_5,'Smoother','off','on',...
        'coarse_solver_flag=SMOOTHER;');
    menu_item(f_mg_5,'PCG','off','off',...
        'coarse_solver_flag = C_PCG;');
    menu_item(f_mg_5,'BiCG-STAB','off','off',...
        'coarse_solver_flag = C_BICG_STAB;');
f_mg_51=menu_item(f_mg_5,'GMRES(k)','off','off',...
    'coarse_solver_flag = C_GMRES;');
f_mg_511=menu_header(f_mg_51,'k = ','on','on');
    menu_item(f_mg_511,'1','off','on','restart=1;');
    menu_item(f_mg_511,'5','off','on','restart=5;');
    menu_item(f_mg_511,'10','off','on','restart=10;');
    menu_item(f_mg_511,'15','off','on','restart=15;');
    menu_item(f_mg_511,'20','off','on','restart=20;');

f_mg_6=menu_header(f_solver_mg,'Coarse-grid Operator','off','on');
    menu_item(f_mg_6,'Standard 5pt','off','on',...
        'coarsening_flag=STANDARD;');
    menu_item(f_mg_6,'Galerkin coarsening','off','off',...
        'coarsening_flag=GALERKIN;');
    menu_item(f_mg_6,'Coeff. Averaging','off','off',...
        'coarsening_flag = AVERAGING;');

f_mg_7=menu_header(f_solver_mg,'MG Cycle','off','on');
    menu_item(f_mg_7,'V-Cycle','off','on',...
        'cycle_flag=V_CYCLE;');
    menu_item(f_mg_7,'W-Cycle','off','on',...
        'cycle_flag=W_CYCLE;');
    menu_item(f_mg_7,'Half V-Cycle','off','off',...
        'cycle_flag=HALF_V_CYCLE;');

% == Results Menu Item ===========================================

f_results=menu_header(main_fig,'Visualize','on','on');

    menu_item(f_results,'Convergence History','off','on',...
        ' subplot(1,1,1);semilogy(its1,resids1,''r-'',its1,resids1,''wo'')');
    menu_item(f_results,'Computed Solution (surf)','off','on',...
        ' subplot(1,1,1);surf(reshape(sol1,nx1,ny1));shading interp;');
    menu_item(f_results,'Computed Solution (pcolor)','off','on',...
        ' subplot(1,1,1);pcolor(reshape(sol1,nx1,ny1));shading interp;');
f_results_1=menu_header(f_results,'X-Axis','off','on');
    menu_item(f_results_1,'Iterations','off','on',...
        'x_axis_flag=ITERATIONS;');
    menu_item(f_results_1,'Time','off','off',...
        'x_axis_flag=TIME;');
    menu_item(f_results_1,'MFlops','off','off',...
        'x_axis_flag=MFLOPS;');
f_results_2=menu_header(f_results,'Y-Axis','off','on');
    menu_item(f_results_2,'Residual','off','on',...
        'y_axis_flag=ITERATIONS;');
    menu_item(f_results_2,'Precon. Residual','off','off',...
        'y_axis_flag=RESIDUAL;');
    menu_item(f_results_2,'MFlops','off','off',...
        'y_axis_flag=PRECON_RESIDUAL;');

f_demos=menu_header(main_fig,'Demos','on','on');
    menu_item(f_demos,'Smoothers','off','on','demo1;');
    menu_item(f_demos,'Fourier analysis','off','on','demo2;');
    menu_item(f_demos,'Truncation error','off','on','demo3;');
