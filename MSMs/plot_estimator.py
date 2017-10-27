from msmbuilder.io import load_generic
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np


def plot_eigenvectors(model, which='left', number=2, fname='msm_eigenvectors'):
    n_states = model.n_states_
    fig, ax = plt.subplots(ncols=1, nrows=number)
    for n_ev in range(number):
        if which == 'left':
            ax[n_ev].bar(range(n_states), model.left_eigenvectors_[:, n_ev])
        else:
            ax[n_ev].bar(range(n_states), model.right_eigenvectors_[:,n_ev])
        ax[n_ev].set_ylabel('{}\n eigenvector'.format(which))
    ax[n_ev].set_xlabel('State')
    plt.savefig('figures/{1}-{0}.png'.format(fname, which))

def plot_single_var(results_df, which, fname='results'):
    """
    plots single variable against the score. 
    :param results_df: appopriately subsetted dataframe
    :param which: which parameter to plot
    :return: None
    """
    label = which.split('__')[1]
    fig, ax = plt.subplots()
    x = results_df[which].values
    y = results_df['mean_test_score']
    err = results_df['std_test_score']
    ax.errorbar(x, y, err)
    ax.set_xscale("log")
    ax.set_ylabel('Score')
    ax.set_xlabel('{}'.format(label))
    plt.savefig('figures/{0}-{1}.png'.format(fname, label))


search_params = load_generic('models/rmsd_model.pickl')
df = pd.DataFrame(search_params.cv_results_)
plot_single_var(df,which='param_cluster__n_clusters' )
print(df.head())
# best_model = search_params.best_estimator_
# msm = best_model.named_steps['msm']
#
# plot_eigenvectors(model=msm, number=3)
# plot_eigenvectors(model=msm, number=3, which='right')



