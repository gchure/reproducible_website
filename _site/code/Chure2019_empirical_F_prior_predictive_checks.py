# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd

n_draws = 1000
n_data = 10
fc_mu = np.random.uniform(0, 1, n_draws)
fc_sigma = np.abs(np.random.normal(0, 0.1, n_draws))

dfs = []
for i in range(n_draws):
    df = pd.DataFrame([] )
    fc_draw = np.random.normal(fc_mu[i], fc_sigma[i], n_data)
    df['fold_change'] = fc_draw
    df['mu'] = fc_mu[i]
    df['sigma'] = fc_sigma[i]
    df['draw'] = i
    df['empirical_F'] = -np.log(fc_mu[i]**-1 - 1)
    dfs.append(df)
df = pd.concat(dfs)
df.to_csv('../../data/Chure2019_empirical_F_prior_predictive_checks.csv', index=False)

