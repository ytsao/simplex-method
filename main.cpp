#include "include/SOPMatrix.hpp"
#include "include/SOPModel.hpp"
#include "include/SOPlpsolver.hpp"

int main() {

#pragma region Mathematical Programming Modeling
    SOPModel model;
    SOPExpression expr;

    // create decision variables
    SOPVariables x1 = model.create_nonbasic("x1", 'F', 0, 0, INT_MAX);
    SOPVariables x2 = model.create_nonbasic("x2", 'F', 1, 0, INT_MAX);

    std::cout << x1.get_var_name() << std::endl;
    std::cout << x2.get_var_name() << std::endl;

    // create objective function
    SOPObjective obj;
    obj.addTerm(-2, x1);
    obj.addTerm(-1, x2);
    model.addMax(obj);

    std::cout << "header: " << model.get_header() << std::endl;

    // create bunch of constraints
    expr.addTerm(-3, x1);
    expr.addTerm(-1, x2);
    model.addLe(expr, -3, "cons1");

    expr.addTerm(-4, x1);
    expr.addTerm(-3, x2);
    model.addLe(expr, -6, "cons2");

    expr.addTerm(-1, x1);
    expr.addTerm(-2, x2);
    model.addLe(expr, -3, "cons3");

    model.create_basis();
    std::cout << "basis: " << model.get_basis() << std::endl;
    std::cout << "all ranges: " << model.get_all_range() << std::endl;
    std::cout << "N: " << model.get_N() << std::endl;
    std::cout << "rhs: " << model.get_rhs() << std::endl;

    for (int i = 0; i < 2; i++) std::cout << "basic variable index: " << model.get_basic(i).get_var_id() << std::endl;

#pragma endregion

#pragma region solve linear program
    SOPlpSolve lpsolver;
    lpsolver.solve(model, 2);
    std::cout << "after solving, rhs: " << model.get_rhs() << std::endl;
    for (int i = 0; i < 2; i++) std::cout << "basic variable index: " << model.get_basic(i).get_var_id() << std::endl;
#pragma endregion

    std::cout << "finish" << std::endl;
    return 0;
}
