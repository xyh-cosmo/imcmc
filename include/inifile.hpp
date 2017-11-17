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

	bool ReadIni( std::string& inifile, bool& save_a_record=true );

	bool ReadInt( std::string& varname, int& value );
	bool ReadBool( std::string& varname, bool& value );
	bool ReadStr( std::string& varname, std::string& value );
	bool ReadDouble( std::string& varname, double& value );

	bool ReadBoolArray( std::string& varname, bool* values, int& nvals );
	bool ReadIntArray( std::string& varname, int* values, int& nvals );
	bool ReadStrArray( std::string& varname, std::string* values, int& nvals );
	bool ReadDoubleArray( std::string& varname, double* values, int& nvals );
};



#endif  //  __IMCMC_INIFILE__
