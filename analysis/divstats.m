%% PFA Divergence plots and statistics

%% MEF MycRas

for tau = 20:5:30
    for bins = [3 5 10 15 20 30]
        myc_exp = ['MEF_MycRas', '_t', num2str(tau), '_b', num2str(bins)];
        myc_exp_name = ['MEF MycRas', ' ', 'tau', num2str(tau), ', ', 'bins', num2str(bins)];
        filename = ['MEF_MycRas_', num2str(tau), '_b', num2str(bins), '_state_vectors.csv'];
        myc_state_vectors = csvread(fullfile('~/src/hm_analysis/data/mef/', filename), 1);
        bin_size = [bins bins];
        pfa_divergence(myc_exp, myc_exp_name, myc_state_vectors, bin_size);
    end
end
%% MEF WT

for tau = 20:5:30
    for bins = [3 5 10 15 20 30]
        wt_exp = ['MEF_WT', '_t', num2str(tau), '_b', num2str(bins)];
        wt_exp_name = ['MEF WT', ' ', 'tau', num2str(tau), ', ', 'bins', num2str(bins)];
        filename = ['MEF_WT_', num2str(tau), '_b', num2str(bins), '_state_vectors.csv'];
        wt_state_vectors = csvread(fullfile('~/src/hm_analysis/data/mef/', filename), 1);
        bin_size = [bins bins];
        pfa_divergence(wt_exp, wt_exp_name, wt_state_vectors, bin_size);
    end
end


%% MuSC FGF2+

for tau = 20:5:30
    for bins = [3 5 10 15 20 30]
        fgf2_exp = ['MuSC_FGF2', '_t', num2str(tau), '_b', num2str(bins)];
        fgf2_exp_name = ['MuSC FGF2', ' ', 'tau', num2str(tau), ', ', 'bins', num2str(bins)];
        filename = ['MuSC_fgf2_', num2str(20), '_b', num2str(bins), '_state_vectors.csv'];
        fgf2_state_vectors = csvread(fullfile('~/src/hm_analysis/data/musc/', filename), 1);
        bin_size = [bins bins];
        pfa_divergence(fgf2_exp, fgf2_exp_name, fgf2_state_vectors, bin_size);
    end
end

%% MuSC FGF-

for tau = 20:5:30
    for bins = [3 5 10 15 20 30]
        nofgf2_exp = ['MuSC_noFGF2', '_t', num2str(tau), '_b', num2str(bins)];
        nofgf2_exp_name = ['MuSC noFGF2', ' ', 'tau', num2str(tau), ', ', 'bins', num2str(bins)];
        filename = ['MuSC_nofgf2_', num2str(20), '_b', num2str(bins), '_state_vectors.csv'];
        nofgf2_state_vectors = csvread(fullfile('~/src/hm_analysis/data/musc/', filename), 1);
        bin_size = [bins bins];
        pfa_divergence(nofgf2_exp, nofgf2_exp_name, nofgf2_state_vectors, bin_size);
    end
end

%% Myoblast FGF2+
for tau = 20:5:30
    for bins = [3 5 10 15 20 30]
        fgf2_exp = ['Myoblast_FGF2', '_t', num2str(tau), '_b', num2str(bins)];
        fgf2_exp_name = ['Myoblast FGF2', ' ', 'tau', num2str(tau), ', ', 'bins', num2str(bins)];
        filename = ['Myoblast_fgf2_', num2str(20), '_b', num2str(bins), '_state_vectors.csv'];
        fgf2_state_vectors = csvread(fullfile('~/src/hm_analysis/data/myoblast/', filename), 1);
        bin_size = [bins bins];
        pfa_divergence(fgf2_exp, fgf2_exp_name, fgf2_state_vectors, bin_size);
    end
end

%% Myoblast FGF2-

for tau = 20:5:30
    for bins = [3 5 10 15 20 30]
        nofgf2_exp = ['Myoblast_noFGF2', '_t', num2str(tau), '_b', num2str(bins)];
        nofgf2_exp_name = ['Myoblast noFGF2', ' ', 'tau', num2str(tau), ', ', 'bins', num2str(bins)];
        filename = ['Myoblast_nofgf2_', num2str(20), '_b', num2str(bins), '_state_vectors.csv'];
        nofgf2_state_vectors = csvread(fullfile('~/src/hm_analysis/data/myoblast/', filename), 1);
        bin_size = [bins bins];
        pfa_divergence(nofgf2_exp, nofgf2_exp_name, nofgf2_state_vectors, bin_size);
    end
end
%% Power Flier to fractal Brownian motion

pwr2fbm_exp = 'pwr2fbm';
pwr2fbm_exp_name = 'Flier to fBm';
pwr2fbm_state_vectors = csvread('~/src/hm_analysis/data/sims/pwr2fbm_state_vectors.csv', 1);
pfa_divergence(pwr2fbm_exp, pwr2fbm_exp_name, pwr2fbm_state_vectors, bin_size);

%% Power Flier to Random Walk

for bins = [3 5 10 15 20 30]
    pwr2rw_exp = ['pwr2rw', '_b', num2str(bins)];
    pwr2rw_exp_name = ['Flier to Random Walk', ' ', 'bins', num2str(bins)];
    filename = ['pwr2rw', '_b', num2str(bins), '_state_vectors.csv'];
    bin_size = [bins bins];
    pwr2rw_state_vectors = csvread(fullfile('~/src/hm_analysis/data/sims/', filename), 1);
    pfa_divergence(pwr2rw_exp, pwr2rw_exp_name, pwr2rw_state_vectors, bin_size);
end 
%% fractal Brownian Motion to Random Walk 

fbm2rw_exp = 'fbm2rw';
fbm2rw_exp_name = 'fBm to Random Walk';
fbm2rw_state_vectors = csvread('~/src/hm_analysis/data/sims/fbm2rw_state_vectors.csv', 1);
pfa_divergence(fbm2rw_exp, fbm2rw_exp_name, fbm2rw_state_vectors, bin_size);

%% Power Flier to Power Flier

pwr2pwr_exp = 'pwr2pwr';
pwr2pwr_exp_name = 'Invariant Flier';
pwr2pwr_state_vectors = csvread('~/src/hm_analysis/data/sims/pwr2pwr_state_vectors.csv', 1);
pfa_divergence(pwr2pwr_exp, pwr2pwr_exp_name, pwr2pwr_state_vectors, bin_size);

%% Random Walk to Random Walk 

for bins = [3 5 10 15 20 30]
    rw2rw_exp = ['rw2rw', '_b', num2str(bins)];
    rw2rw_exp_name = ['Invariant Random Walk', ' ', 'bins', num2str(bins)];
    filename = ['rw2rw', '_b', num2str(bins), '_state_vectors.csv'];
    rw2rw_state_vectors = csvread(fullfile('~/src/hm_analysis/data/sims/', filename), 1);
    pfa_divergence(rw2rw_exp, rw2rw_exp_name, rw2rw_state_vectors, bin_size);
end