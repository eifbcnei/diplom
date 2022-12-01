run('set_up_global_config.m')
disp('Fourier')
run('fourier_scan.m')
logger = '';
disp('Capon')
run('capon_scan.m')
logger = '';
disp('Heat noise')
run('heat_noise_scan.m')
%run('global_sandbox.m')
run('print_results.m')