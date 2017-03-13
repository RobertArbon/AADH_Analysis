# Compares the results of manual CV grid search and sklearn CV search for understanding

import pandas as pd

manual = pd.read_pickle('manual_results.pickl')
sklrn = pd.read_pickle('sklrn_results.pickl')

manual['err'] *= 0.5

manual.rename(columns={'n_componets': 'n_components',
                       'err': 'std'}, inplace=True)
sklrn.rename(columns={'mean_test_score': 'mean',
                      'std_test_score': 'std',
                      'param_lag_time': 'lag_time',
                      'param_n_components': 'n_components'}, inplace=True)

sklrn = sklrn[['mean', 'std', 'lag_time', 'n_components']]

df = pd.merge(left=manual, right=sklrn, on=['lag_time', 'n_components'], how='inner')

df['mean_diff'] = 100*(df['mean_x']/df['mean_y']-1)
df['std_diff'] = 100*(df['std_x']/df['std_y']-1)

print(df)


# This ran in fit_tica_grid_search.py to check whether the CV strategy was
# running as I thought it was.
#
# results = {'n_componets': [],
#            'lag_time': [],
#            'mean': [],
#            'err': []}
#
# for comp in params['n_components']:
#     for lag in params['lag_time']:
#         temp = []
#         for train_idx, test_idx in cv_iter.split(X):
#             model = tICA(kinetic_mapping=True, lag_time=lag, n_components=comp)
#             train_set = [X[i] for i in train_idx]
#             test_set = [X[i] for i in test_idx]
#             model.fit(train_set)
#             temp.append(model.score(test_set))
#         mean = np.array(temp).mean()
#         std = np.array(temp).std()
#         results['n_componets'].append(comp)
#         results['lag_time'].append(lag)
#         results['mean'].append(mean)
#         results['err'].append(std * 2)
#
# df = pd.DataFrame(results)
# save_generic(df, 'manual_results.pickl')
# print(df)