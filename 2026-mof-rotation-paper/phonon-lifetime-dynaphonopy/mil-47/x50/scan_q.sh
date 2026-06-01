
dynaphopy input_file mil-47.lammpstrj -sp mode_spectrum_0 -q 0 0 0 -r 0 100 -psm 2 -pa --silent > fitting_0.log

fitdata mode_spectrum_0 --silent > fitting_0.log

dynaphopy input_file mil-47.lammpstrj -sp mode_spectrum_1 -q 0.125 0 0 -r 0 100 -psm 2 -pa --silent > fitting_1.log

fitdata mode_spectrum_1 --silent > fitting_1.log

dynaphopy input_file mil-47.lammpstrj -sp mode_spectrum_2 -q 0.25 0 0 -r 0 100 -psm 2 -pa --silent > fitting_2.log

fitdata mode_spectrum_2 --silent > fitting_2.log

dynaphopy input_file mil-47.lammpstrj -sp mode_spectrum_3 -q 0.375 0 0 -r 0 100 -psm 2 -pa --silent > fitting_3.log

fitdata mode_spectrum_3 --silent > fitting_3.log

dynaphopy input_file mil-47.lammpstrj -sp mode_spectrum_4 -q 0.5 0 0 -r 0 100 -psm 2 -pa --silent > fitting_4.log

fitdata mode_spectrum_4 --silent > fitting_4.log

dynaphopy input_file mil-47.lammpstrj -sp mode_spectrum_5 -q 0.625 0 0 -r 0 100 -psm 2 -pa --silent > fitting_5.log

fitdata mode_spectrum_5 --silent > fitting_5.log

dynaphopy input_file mil-47.lammpstrj -sp mode_spectrum_6 -q 0.75 0 0 -r 0 100 -psm 2 -pa --silent > fitting_6.log

fitdata mode_spectrum_6 --silent > fitting_6.log

dynaphopy input_file mil-47.lammpstrj -sp mode_spectrum_7 -q 0.875 0 0 -r 0 100 -psm 2 -pa --silent > fitting_7.log

fitdata mode_spectrum_7 --silent > fitting_7.log
