#!/bin/bash -l
'/share/pkg.7/matlab/2020b/install/bin/matlab' -nodesktop -noFigureWindows -nosplash -singleCompThread -logfile '/project/busplab/software/ecog/classification_AM/.qlog/210517002430258/node.0001210517002430258.stdlog' -r "addpath '/project/busplab/software/spm12'; addpath '/project/busplab/software/conn'; cd '/project/busplab/software/ecog/classification_AM/.qlog/210517002430258'; conn_jobmanager('rexec','/project/busplab/software/ecog/classification_AM/.qlog/210517002430258/node.0001210517002430258.mat'); exit"
echo _NODE END_
