import logging

from lpsolve55 import lpsolve
import pyomo.environ as pyo

from rmgpy.exceptions import ILPSolutionError


def solve_clar_lp_pyomo(objective, a, l, m, n, exo, constraints, mol):
    # Setup the LP Problem using pyomo
    lp_model = pyo.ConcreteModel()
    lp_model.i = pyo.RangeSet(0, m-1)
    lp_model.j = pyo.RangeSet(0, n-1)
    lp_model.x = pyo.Var(lp_model.j, domain=pyo.Binary)
    lp_model.a = pyo.Param(lp_model.i, lp_model.j, initialize=lambda _, i, j: a[i][j])
    lp_model.o = pyo.Param(lp_model.j, domain=pyo.Binary, initialize=lambda _, j: objective[j])

    def objective_function(model):
        return pyo.summation(model.x, model.o, index=model.j)

    lp_model.obj = pyo.Objective(rule=objective_function, sense=pyo.maximize)

    def constraint_function(model, i):
        return sum(model.x[j]*model.a[i, j] for j in model.j) == 1

    lp_model.equality_constraints = pyo.Constraint(lp_model.i, rule=constraint_function)

    # Constrain values of exocyclic bonds, since we don't want to modify them
    for k in range(l, n):
        if exo[k - l] is not None:
            lp_model.x[k].fix(exo[k - l])

    # Add constraints to problem if provided
    if constraints is not None:
        lp_model.leq_constraints = pyo.ConstraintList()
        for constraint in constraints:
            try:
                lp_model.leq_constraints.add(sum(lp_model.x[j] * constraint[0][j] for j in lp_model.j) <= constraint[1])
            except Exception as e:
                logging.debug('Unable to add constraint: {0} <= {1}'.format(constraint[0], constraint[1]))
                logging.debug(mol.toAdjacencyList())
                if str(e) == 'invalid vector.':
                    raise ILPSolutionError('Unable to add constraint, likely due to '
                                           'inconsistent aromatic ring perception.')
                else:
                    raise e

    opt = pyo.SolverFactory('glpk')
    result = opt.solve(lp_model)

    status = result.solver.status == pyo.SolverStatus.ok
    obj_val = pyo.value(lp_model.obj)
    solution = lp_model.x.extract_values().values()

    return status, obj_val, solution


def solve_clar_lp_lpsolve(objective, a, l, m, n, exo, constraints, mol):
    # Solve LP problem using lpsolve
    lp = lpsolve('make_lp', m, n)  # initialize lp with constraint matrix with m rows and n columns
    lpsolve('set_verbose', lp, 2)  # reduce messages from lpsolve
    lpsolve('set_obj_fn', lp, objective)  # set objective function
    lpsolve('set_maxim', lp)  # set solver to maximize objective
    lpsolve('set_mat', lp, a)  # set left hand side to constraint matrix
    lpsolve('set_rh_vec', lp, [1] * m)  # set right hand side to 1 for all constraints
    lpsolve('set_constr_type', lp, ['='] * m)  # set all constraints as equality constraints
    lpsolve('set_binary', lp, [True] * n)  # set all variables to be binary

    # Constrain values of exocyclic bonds, since we don't want to modify them
    for i in range(l, n):
        if exo[i - l] is not None:
            # NOTE: lpsolve indexes from 1, so the variable we're changing should be i + 1
            lpsolve('set_bounds', lp, i + 1, exo[i - l], exo[i - l])

    # Add constraints to problem if provided
    if constraints is not None:
        for constraint in constraints:
            try:
                lpsolve('add_constraint', lp, constraint[0], '<=', constraint[1])
            except Exception as e:
                logging.debug('Unable to add constraint: {0} <= {1}'.format(constraint[0], constraint[1]))
                logging.debug(mol.toAdjacencyList())
                if str(e) == 'invalid vector.':
                    raise ILPSolutionError('Unable to add constraint, likely due to '
                                           'inconsistent aromatic ring perception.')
                else:
                    raise e

    status = lpsolve('solve', lp)
    if status:
        status = 0

    obj_val, solution = lpsolve('get_solution', lp)[0:2]
    lpsolve('delete_lp', lp)  # Delete the LP problem to clear up memory

    return status, obj_val, solution


solve_clar_lp = solve_clar_lp_lpsolve
# solve_clar_lp = solve_clar_lp_pyomo
