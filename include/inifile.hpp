/*
 *  New parser class for imcmc.
 *  -- Can be viewed as a clone of inifile included in cosmomc which is written in fortran.
 */

#ifndef __IMCMC_INIFILE__
#define __IMCMC_INIFILE__

#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <utility>
#include <map>
#include <cmath>
#include <ctime>

#include <string>
#include <vector>
#include <cstdlib>
#include <stdexcept>

struct ArrayInt {
    std::string         VarName;
    int                 VarNum;
    std::vector<int>    Value;
    int getValue(int idx=0);
};

struct ArrayBool {
    std::string         VarName;
    int                 VarNum;
    std::vector<bool>   Value;
    bool getValue(int idx=0);
};

struct ArrayDouble {
    std::string         VarName;
    int                 VarNum;
    std::vector<double> Value;
    double getValue(int idx=0);
};

struct ArrayString {
    std::string         VarName;
    int                 VarNum;
    std::vector<std::string>    Value;
    std::string getValue(int idx=0);
};

class IniFile {
  private:
    std::vector<bool> Bool;
    std::vector<int> Int;
    std::vector<double> Double;
    std::vector<std::string> String;

    std::vector<std::vector<bool>> BoolArray;
    std::vector<std::vector<int>> IntArray;
    std::vector<std::vector<double>> DoubleArray;
    std::vector<std::vector<std::string>> StringArray;

  public:
    bool    GetBool(std::string pname);
    int     GetInt(std::string pname);
    double  GetDouble(std::string pname);
    std::string  GetString(std::string pname);

    void    GetValue(std::string pname, bool& value);
    void    GetValue(std::string pname, int& value);
    void    GetValue(std::string pname, double& value);
    void    GetValue(std::string pname, std::string& value);
    void    GetValue(std::string pname, bool *value, int n_value);
    void    GetValue(std::string pname, int *value, int n_value);
    void    GetValue(std::string pname, double *value, int n_value);
    void    GetValue(std::string pname, )
};



#endif  //  __IMCMC_INIFILE__
