"""Independent SciPy/HiGHS check of the R/lpSolve guarded component-CVaR LP."""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.optimize import linprog
from scipy import sparse

root = Path('outputs')
d = pd.read_csv(root / 'LP_CROSSCHECK_input.csv')
q = d.prob.to_numpy()
returns = d.drop(columns='prob').to_numpy()
m, p = returns.shape
size = p + 1 + 2 + 2*m + 1 + m
objective = np.zeros(size)
objective[p] = 1
rows, rhs = [], []
for k, scale in enumerate([1, 1.5]):
    z, u = p + 1 + k, p + 3 + k*m
    for j in range(m):
        row = np.zeros(size)
        row[:p], row[z], row[u+j] = -returns[j]*scale, -1, -1
        rows.append(row)
        rhs.append(0)
    row = np.zeros(size)
    row[z], row[u:u+m], row[p] = 1, q/.05, -1
    rows.append(row)
    rhs.append(0)
z = p + 3 + 2*m
u = z + 1
for j in range(m):
    row = np.zeros(size)
    row[:p], row[z], row[u+j] = -returns[j], -1, -1
    rows.append(row)
    rhs.append(0)
expected = pd.read_csv(root / 'LP_CROSSCHECK_expected.csv')
row = np.zeros(size)
row[z], row[u:u+m] = 1, q/.05
rows.append(row)
rhs.append(expected.guard.iloc[0])
bounds = [(0,.4)]*p + [(None,None)]*3 + [(0,None)]*(2*m) + [(None,None)] + [(0,None)]*m
eq = np.zeros((1,size))
eq[0,:p] = 1
fit = linprog(objective, A_ub=sparse.csr_matrix(rows), b_ub=rhs,
              A_eq=eq, b_eq=[1], bounds=bounds, method='highs')
assert fit.success, fit.message
objective_error = abs(fit.fun-expected.risk.iloc[0])
weight_error = np.max(abs(fit.x[:p]-expected.weight))
assert objective_error < 1e-8 and weight_error < 1e-6
pd.DataFrame([dict(scipy_success=fit.success,objective_error=objective_error,
                  weight_error=weight_error,pass_check=True)]).to_csv(root/'INDEPENDENT_LP_CHECK.csv',index=False)
print('Independent HiGHS check passed',objective_error,weight_error)
