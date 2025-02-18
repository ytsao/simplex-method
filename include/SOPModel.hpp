#ifndef _SOPMODEL_H_
#define _SOPMODEL_H_

#include "SOPMatrix.hpp"

class SOPModel {
private:
    int num_rows = 0;
    int num_cols = 0;
    std::vector<SOPVariables> nonbasic;
    std::vector<SOPVariables> basic;
    SOPObjective obj;
    std::vector<SOPRange> cons;
    SOPVector b;
    SOPVector c;
    SOPMatrix A;
    SOPMatrix N;
    SOPMatrix B;
    SOPMatrix eta;

protected:
public:
    SOPModel() {
    }

    SOPVariables create_nonbasic(std::string name, char type, int id, double lb, double ub){
      SOPVariables new_var = SOPVariables(name, type, id, lb, ub);
      this->nonbasic.push_back(new_var);
      this->num_cols;
      
      return new_var;
    }

    void addMax(SOPObjective _obj){
      this->c = _obj.get_all_coef();
      for (int i = 0; i < _obj.getSize(); ++i){
        int tmp = this->c.coeff(i);
        this->c(i) = -tmp;
      }
      _obj.clear();
    }

    void addMin(SOPObjective _obj){
      this->c = _obj.get_all_coef();
      _obj.clear();
    }

    SOPRange addLe(SOPExpression& con, double rhs, std::string con_name){
      this->A.conservativeResize(this->A.rows()+1, get_num_nonbasic());
      this->N.conservativeResize(this->N.rows()+1, get_num_nonbasic());

      for (int i = 0; i < get_num_nonbasic(); ++i){
        this->A(this->A.rows()-1, con.get_var(i).get_var_id()) = con.get_coef(i);
        this->N(this->N.rows()-1, con.get_var(i).get_var_id()) = con.get_coef(i);
      }

      this->b.conservativeResize(this->b.size()+1);
      this->b(this->b.size()-1) = rhs;

      // create basic variable
      SOPVariables basic = SOPVariables("x" + std::to_string(this->get_num_cols()+this->get_num_rows()+1), 'R', this->get_num_cols()+this->get_num_rows(), 0, INT_MAX);
      this->basic.push_back(basic);

      // create range
      SOPRange new_con = SOPRange(con.get_all_coef(), con.get_all_vars(), rhs, 'L', con_name);
      this->cons.push_back(new_con);
      this->num_rows++;
      con.clear();

      return new_con;
    }

    SOPRange addGe(SOPExpression& con, double rhs, std::string con_name){
      this->A.conservativeResize(this->A.rows()+1, get_num_nonbasic());
      this->N.conservativeResize(this->N.rows()+1, get_num_nonbasic());

      for (int i = 0; i < get_num_nonbasic(); ++i){
        this->A(this->A.rows()-1, con.get_var(i).get_var_id()) = -con.get_coef(i);
        this->N(this->N.rows()-1, con.get_var(i).get_var_id()) = -con.get_coef(i);
      }

      this->b.conservativeResize(this->b.size()+1);
      this->b(this->b.size()-1) = rhs;

      // create basic variable
      SOPVariables basic = SOPVariables("x"+std::to_string(this->get_num_cols()+this->get_num_rows()+1), 'R', this->get_num_cols()+this->get_num_rows(), 0, INT_MAX);
      this->basic.push_back(basic);

      // create range
      SOPRange new_con = SOPRange(con.get_all_coef(), con.get_all_vars(), rhs, 'G', con_name);
      this->cons.push_back(new_con);
      this->num_rows++;
      con.clear();

      return new_con;
    }

    SOPRange addEq(SOPExpression& con, double rhs, std::string con_name){
      this->A.conservativeResize(this->A.rows()+1, get_num_nonbasic());
      this->N.conservativeResize(this->N.rows()+1, get_num_nonbasic());

      for (int i = 0; i < get_num_nonbasic(); ++i){
        this->A(this->A.rows()-1, con.get_var(i).get_var_id()) = con.get_coef(i);
        this->N(this->N.rows()-1, con.get_var(i).get_var_id()) = con.get_coef(i);
      }

      this->b.conservativeResize(this->b.size()+1);
      this->b(this->b.size()-1) = rhs;

      // create range
      SOPRange new_con = SOPRange(con.get_all_coef(), con.get_all_vars(), rhs, 'E', con_name);
      this->cons.push_back(new_con);
      this->num_rows++;
      con.clear();

      return new_con;
    }

    void create_basis(){
      this->B = SOPMatrix::Identity(this->num_rows, this->num_rows);
      this->A.conservativeResize(this->A.rows(), this->nonbasic.size()+this->basic.size());
      this->A.leftCols(this->N.cols()) = this->N;
      this->A.rightCols(this->B.cols()) = this->B;
    }

    int get_num_rows(){
      return this->num_rows;
    }

    int get_num_cols(){
      return this->num_cols;
    }

    int get_num_nonbasic(){
      return this->nonbasic.size();
    }

    int get_num_basic(){
      return this->basic.size();
    }

    int get_all_variables(){
      return this->nonbasic.size()+this->basic.size();
    }

    SOPVariables get_nonbasic(int position){
      return this->nonbasic[position];
    }

    SOPVariables get_basic(int position){
      return this->basic[position];
    }

    SOPVector get_header(){
      return this->c;
    }

    SOPVector get_rhs(){
      return this->b;
    }

    SOPMatrix get_all_range(){
      return this->A;
    }

    SOPMatrix get_basis(){
      return this->B;
    }

    SOPMatrix get_N(){
      return this->N;
    }

    void set_header(int position, double value){
      this->c(position) = value;
    }

    void set_rhs(int position, double value){
      this->b(position) = value;
    }

    void set_nonbasic(int position, SOPVariables var){
      this->nonbasic[position] = var;
    }

    void set_basic(int position, SOPVariables var){
      this->basic[position] = var;
    }

    void set_eta(int leaving_index, std::vector<double> delta_x_B){
      this->eta.resize(this->num_rows, this->num_rows);
      this->eta.setZero();
      this->eta.diagonal().array() += 1;
      for (int i = 0; i < this->num_rows; ++i){
        double value = 0;
        if (i != leaving_index){
          value = -(delta_x_B[i] / delta_x_B[leaving_index]);
        }
        else{
          value = 1 / delta_x_B[leaving_index];
        }
        this->eta(i, leaving_index) = value;
      }
      std::cout << "create identical matrix" << std::endl;
    }

    void update_basis(){
      this->B = this->eta * this->B;
    }

    void clear();
};

class SOPSolver {
private:
protected:
public:
    SOPSolver() {
    }

    bool solve(SOPModel model, int algorithm){
      return true;
    }
};

#endif
