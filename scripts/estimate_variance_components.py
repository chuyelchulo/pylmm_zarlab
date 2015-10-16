import numpy as np
import numpy.linalg as npl
import itertools
from pylmm import input, lmm


def simplex(dimensions, fineness):
    possibles = itertools.product(*[np.linspace(0, 1, fineness+1)]*dimensions)
    goods = (x for x in possibles if np.isclose(sum(x), 1))
    return goods


def log_likelihood(residuals, V, V_inv):
    n_indivs = len(Y)
    indiv_part = -.5 * n_indivs * np.log(2*np.pi)
    det_part = -.5 * np.log(npl.det(V))
    residual_part = -.5 * np.dot((residuals.T, V_inv), residuals)
    LL = indiv_part + det_part + residual_part
    return LL


def grid_search(K_G, K_GxE, X, Y, fineness):
    assert K_G.shape == K_GxE.shape
    I = np.eye(K_G.shape[0])
    best_LL = None
    best_var_comps = None

    for var_comps in simplex(3, fineness):
        Var_G, Var_GxE, Var_E = var_comps
        V = Var_G*K_G + Var_GxE * K_GxE + Var_E * I
        V_inv = npl.inv(V)
        X_pinv = npl.inv(np.dot(np.dot(X.T, V_inv), X))
        rotated_dot_product = np.dot(np.dot(X.T, V_inv), Y)
        betas = np.dot(X_pinv, rotated_dot_product)
        residuals = Y - np.dot(X, betas)
        LL = log_likelihood(residuals, V, V_inv)

        if np.isnan(LL):
            continue

        if best_LL is None:
            best_LL = LL
            best_var_comps = var_comps
        else:
            if LL > best_LL:
                best_LL = LL
                best_var_comps = var_comps
            else:
                pass
    return best_var_comps


def impute_columns(X):
    for i in range(X.shape[1]):
        col = X[:, i]
        mean_val = np.nanmean(col)
        X[np.isnan(col), i] = mean_val
    return X


def make_K_GxE(K_G, env):
    GxE_mask = np.array([[1 if env1 == env2 else 0 for env2 in env]
                         for env1 in env])
    K_GxE = K_G*GxE_mask
    return K_GxE


def main(Y, X, K_G, env, fineness):
    good_Ys = ~np.isnan(Y)
    good_envs = ~np.isnan(env)
    indiv_mask = good_Ys & good_envs

    good_Y = Y[indiv_mask]
    good_X = X[indiv_mask, :]

    imputed_X = impute_columns(good_X)
    full_X = np.column_stack((imputed_X, env))
    K_GxE = make_K_GxE(K_G, env)
    components = grid_search(K_G=K_G, K_GxE=K_GxE, X=full_X, Y=good_Y, fineness=fineness)
    return components


if __name__ == "__main__":
    # v = simplex(3, 10)
    # for u in v:
    #     print u

    X = np.eye(6)
    X[1:4, 2:6] = np.nan
    print X
    print impute_columns(X)
