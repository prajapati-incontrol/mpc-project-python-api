clc, clear

x = input('1) Run Full-state feedback MPC \n 2) Run output feedback disturbance rejection MPC \n 3) Run the stability analysis: Plot CLF decrease \n 4) Plot the X_N estimate\n');
switch x 
    case 1
        clc, clear, close all
        run('regulator_MPC_V9.m')
    case 2
        clc, clear, close all
        run('Output_Feedback_Disturbance_MPC.m')
    case 3 
        clc, clear, close all
        run('stability_analysis.m')
    case 4 
        clc, clear, close all
        % Note: this code takes about 82 minutes!
        run("Final_Xn_estimation.m")
    otherwise
        disp('\n\n Please select a valid input \n\n Run the main.m again');
end
