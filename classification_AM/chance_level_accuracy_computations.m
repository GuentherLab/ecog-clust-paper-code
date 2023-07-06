%     calculate chance-level decoding accuracy for consonants, vowels, and syllables for each subject
%
%%% script adapted from code supplemental to paper (combrisson and jerbi 2015): 
%%%       'Exceeding chance level by chance: The caveat of theoretical chance levels in brain signal classification and statistical assessment of decoding accuracy'

%% 
% Description:
% This code computes and plots the statistical chance level as a function of trial number using the cumulative
% binomial distribution function.
%
%
% by:
% Etienne Combrisson (1,2) [PhD student] / Contact: etienne.combrisson@inserm.fr 
% Karim Jerbi (1,3) [PhD, Assistant Professor] 
% 1 DYCOG Lab, Lyon Neuroscience Research Center, INSERM U1028, UMR 5292, University Lyon I, Lyon, France
% 2 Center of Research and Innovation in Sport, Mental Processes and Motor Performance, University of Lyon I, Lyon, France
% 3 Psychology Department, University of Montreal, QC, Canada
%
% Version 1.0
% Created: 21/10/2014
% Latest update: - 23/10/14
%

% adapted 2023/6/30 by andrew meier for paper: "Lateralization and time-course of cortical phonological representations during syllable production" 

clear
set_paths()


%% params

pvals_to_test = 0.05; %%% can be vector of multiple values
% pvals_to_test = 1e-7; %%% can be vector of multiple values

%% import subject info, organize variables
dat_2sec = load([DIR_LDA_2SEC_WINDOW filesep 'single_elc_included_trialwise_accuracy.mat']); 
subdata = dat_2sec.subdata; % includes individual-trial performance for each subject
% % % % subdata = readtable(SUBJECT_TABLE_FILENAME);
nsubs = height(subdata);
phonunits =              {'cons','vow','syl'};
    nlabels_per_phonunit = [4,     3,   12  ]; 
nphonunits = length(phonunits);

% add variables
nancol = nan(nsubs,1);
subdata = [subdata, table(nancol,nancol,nancol,'VariableNames',{'cons_bino_chance','vow_bino_chance','syl_bino_chance'})];

nancol = nan(nphonunits,1);
celcol = cell(nphonunits,1);
phontab = table(phonunits',nlabels_per_phonunit',nancol,nancol,nancol,celcol,'VariableNames',{'phonunit','nlabels','acc_mean','acc_chance','ntrials','trials'});

%% compute chance level accuracies for each subject
for isub = 1:nsubs
    ntrials = subdata.ntrials(isub); % Number of epoch
    for iphon = 1:nphonunits
        thisphon = phonunits{iphon};
        nlabels = nlabels_per_phonunit(iphon); % Number of classes

        %  Construction of label vectors
        %%%%%%%% make a vector ytemp with a label value for each trial; this vector will be supplied to the function StatTh()
        ytemp = [];
        for k=1:nlabels
            ytemp = [ytemp;k.*ones(round(ntrials/nlabels),1)];
        end
        
        % Construction of stat values:
        for ipval = 1:length(pvals_to_test) % p value loop
            pStr{ipval} = ['p<' num2str(pvals_to_test(ipval))];
            sgnf_thresh(ipval) = StatTh(ytemp,pvals_to_test(ipval)) / 100;
        end

        % use the chance level for the first specified p value
        subdata{isub, [thisphon, '_bino_chance']} = sgnf_thresh(1);
    end
end

%% compute chance level accuracies for all subjects combined
for iphon = 1:nphonunits
    thisphon = phonunits{iphon};
    nlabels = nlabels_per_phonunit(iphon); % Number of classes

    % concatanate accuracies from all trials from all subjects for this phonunit
    phontab.trials{iphon} = cell2mat(cellfun(@transpose,subdata{:,[thisphon,'_acc_all_elcs']},'UniformOutput',false));
    phontab.acc_mean(iphon) = mean(phontab.trials{iphon}); 
    phontab.ntrials(iphon) = length(phontab.trials{iphon}); 

    %  Construction of label vectors
    %%%%%%%% make a vector ytemp with a label value for each trial; this vector will be supplied to the function StatTh()
    ytemp = [];
    for k=1:nlabels
        ytemp = [ytemp;k.*ones( round(phontab.ntrials(iphon)/phontab.nlabels(iphon)) ,1)];
    end
    
    % Construction of stat values:
    for ipval = 1:length(pvals_to_test) % p value loop
        pStr{ipval} = ['p<' num2str(pvals_to_test(ipval))];
        sgnf_thresh(ipval) = StatTh(ytemp,pvals_to_test(ipval)) / 100;
    end

    % use the chance level for the first specified p value
    phontab.acc_chance(iphon) = sgnf_thresh(1);

end







%%
function [threshold,nbclass] = StatTh(y,alpham)
%--------------------------------------------------------------------------
% This function is used to compute a statistical threshold using the
% binomial inverse cumulative distribution
% -> INPUT VARIABLES:
% - y: label vector (Ex: y = [1 1 1 2 2 2 3 3 3 4 4 4];)
% - alpham: p-value (Ex: alpham = 0.01) 
% -> OUTPUT VARIABLES:
% - threshold: as a percentage used for decoding accuracy
% - nbclass: number of class in label vector y
%
% It calls Matlab's built-in "binoinv" function
%
% by:
% Etienne Combrisson (1,2) [PhD student] / Contact: etienne.combrisson@inserm.fr 
% Karim Jerbi (1,3) [PhD, Assistant Professor] 
% 1 DYCOG Lab, Lyon Neuroscience Research Center, INSERM U1028, UMR 5292, University Lyon I, Lyon, France
% 2 Center of Research and Innovation in Sport, Mental Processes and Motor Performance, University of Lyon I, Lyon, France
% 3 Psychology Department, University of Montreal, QC, Canada
%
% Version 1.0
% Created: 21/10/2014
% Latest update: - 23/10/14
%
%--------------------------------------------------------------------------
if (nargin < 2)
    alpham = 0.01;
end
nbepoch = length(y);
nbclass = length(unique(y));
threshold = binoinv(1-alpham,nbepoch,1/nbclass)*100/nbepoch;

end






