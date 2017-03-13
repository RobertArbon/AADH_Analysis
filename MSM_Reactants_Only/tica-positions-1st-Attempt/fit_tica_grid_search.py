#
# Does grid search for optimum parameters
#
from msmbuilder.io import load_trajs, save_trajs, save_generic
from msmbuilder.decomposition import tICA
from sklearn.model_selection import ShuffleSplit, GridSearchCV
import pandas as pd


if __name__ == '__main__':

    # Load data
    feature_name = 'Positions'
    meta, ftrajs = load_trajs("Scaled-{}-ftraj".format(feature_name))
    X = list(ftrajs.values())

    # Specify CV strategy and parameters
    cv_iter = ShuffleSplit(n_splits=10, test_size=0.5, random_state=0)
    param_grid = [{'n_components': [10,20,40],
                   'lag_time': [1,10,100]}]

    # CV object
    model = tICA(kinetic_mapping=True)

    # Do grid search
    clf = GridSearchCV(estimator=model, param_grid=param_grid, cv=cv_iter, n_jobs=2)
    clf.fit(X)

    # Save results
    results = pd.DataFrame(clf.cv_results_)
    save_generic(results, '{}-grid-search-results.pickl'.format(feature_name))

    # Print Results
    print("Best parameters set found on development set:")
    print(clf.best_params_)
    print()
    print("Grid scores on development set:")
    print()
    means = clf.cv_results_['mean_test_score']
    stds = clf.cv_results_['std_test_score']
    for mean, std, params in zip(means, stds, clf.cv_results_['params']):
        print("%0.3f (+/-%0.03f) for %r"
              % (mean, std * 2, params))

    # Fit best estimator to data
    tica = clf.best_estimator_
    ttrajs = {}
    for k, v in ftrajs.items():
        ttrajs[k] = tica.partial_transform(v)

    # Save
    save_trajs(ttrajs, '{}-ttrajs'.format(feature_name), meta)
    save_generic(tica, '{}-tica.pickl'.format(feature_name))

