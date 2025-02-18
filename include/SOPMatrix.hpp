#ifndef _SOPMATRIX_H_
#define _SOPMATRIX_H_

#include <bits/stdc++.h>
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Core>

typedef Eigen::VectorXd SOPVector;
typedef Eigen::MatrixXd SOPMatrix;

class SOPVariables {
private:
    std::string name;
    char type;
    int id;
    double lb;
    double ub;

protected:
public:
    SOPVariables() {
    }

    SOPVariables(std::string _name, char _type, int _id){
      this->name = _name;
      this->type = _type;
      this->id = _id;
      this->lb = 0;
      this->ub = INT_MAX;
    }

    SOPVariables(std::string _name, char _type, int _id, double _lb, double _ub){
      this->name = _name;
      this->type = _type;
      this->id = _id;
      this->lb = lb;
      this->ub = ub;
    }

    std::string get_var_name(){
      return this->name;
    }

    int get_var_id(){
      return this->id;
    }
};

class SOPExpression {
private:
    SOPVector coef;
    std::vector<SOPVariables> vars;

protected:
public:
    SOPExpression() {
    }

    void addTerm(const double coef, const SOPVariables var){
      this->coef.conservativeResize(this->coef.size()+1);
      this->coef(this->coef.size()-1) = coef;
      this->vars.push_back(var);
    }

    void addTerms(const double* coefs, const SOPVariables* vars){
      int size = sizeof(coef) / sizeof(double);
      for (int i = 0; i < size; ++i){
        this->coef.conservativeResize(this->coef.size()+1);
        this->coef(this->coef.size()-1) = coef[i];
        this->vars.push_back(vars[i]);
      }
    }

    int getSize(){
      return this->vars.size();
    }

    double get_coef(const int position){
      return this->coef.coeff(position);
    }

    SOPVector get_all_coef(){
      return this->coef;
    }

    SOPVariables get_var(const int position){
      return this->vars[position];
    }

    std::vector<SOPVariables> get_all_vars(){
      return this->vars;
    }

    void clear(){
      this->coef.resize(0);
      this->vars.clear();
      this->vars.shrink_to_fit();
    }
};

class SOPObjective : public SOPExpression {
private:
    SOPVector coef;
    std::vector<SOPVariables> vars;

protected:
public:
    SOPObjective() {
    }
};

class SOPRange : public SOPExpression {
private:
    SOPVector coef;
    std::vector<SOPVariables> vars;
    double rhs;
    char direction;
    std::string con_name;

protected:
public:
    SOPRange() {
    }

    SOPRange(SOPVector _coef, std::vector<SOPVariables> _vars, double _rhs, char _direction, std::string _con_name){
      this->coef = _coef;
      this->vars = _vars;
      this->rhs = _rhs;
      this->direction = _direction;
      this->con_name = _con_name;
    }
};

#endif
