/*
 *    parser++: a C++ parser
 *    This parser can be embbed into any other projects, it depends only on the standard libs
 *
 *  TODO: some functions of this parser will be re-written in templates.
 */

#ifndef __IMCMC_PARSER__
#define __IMCMC_PARSER__

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


// //  Detector warnings
// #define DetectWarning(condition,msg)                        \
// {                                                           \
//     if( condition == true ) {                               \
//         std::cout << "\n==> Warning detected: \n"           \
//              << "==> File: " << __FILE__ << "\n"            \
//              << "==> Line: " << __LINE__ << "\n"            \
//              << "==> Func: " << __FUNCTION__ << "\n"        \
//              << "==> Warning Info: " << msg << std::endl;   \
//     }                                                       \
// }


// //  Detector error for given condition
// #define DetectError(condition,msg)                      \
// {                                                       \
//     if( condition == true ) {                           \
//         std::cout << "\n==> Error detected: \n"         \
//              << "==> File: " << __FILE__ << "\n"        \
//              << "==> Line: " << __LINE__ << "\n"        \
//              << "==> Func: " << __FUNCTION__ << "\n"    \
//              << "==> Error Info: " << msg << std::endl; \
//     }                                                   \
// }

// //  stop on given condition
// #if defined(__IMCMC_MPI__)

// #define StopOnError(condition,msg)      \
// {                                       \
//     DetectError(condition,msg);         \
//     std::cout << "==> Stop Now!\n";     \
// 	MPI::COMM_WORLD.Abort(MPI::COMM_WORLD.Get_rank());	\
// }

// #else

// #define StopOnError(condition,msg)		\
// {										\
// 	DetectError(condition,msg);			\
// 	std::cout << "==> Stop Now!\n";		\
// 	exit(0);							\
// }

// #endif


namespace imcmc{
namespace parser{

    namespace Read{

        bool Has_File( std::string paramfile );   //  test whether a file exist
        bool Is_Commented( std::string line );    //  check whether the given line is commented, if commented, true will be return

        std::string Remove_Comments( std::string line );  //    extract un-commented parts of a line
        void        Remove_TabSpace( std::string& value, bool info=false );

        bool        Is_Array( std::string param );
        std::string GetArrayName( std::string param );     //    A(:) means that A is a vector, but the size is unknown ...

        //  check if the given line has the key you wanted, not case sensitive
        bool Has_Key( std::string line, std::string key );
        bool Has_Key_in_File( std::string file, std::string key );

        //  find the beginning and ending postions for value(s)
        std::string::size_type Begin_of_Value( std::string line );
        std::string::size_type End_of_Value( std::string line );

        //  count the number of values in the given line
        int Num_of_Value(   std::string            line,
                            std::string::size_type begin,
                            std::string::size_type end,
                            std::string            type="string" );

        int Num_of_Value_for_Key( std::string infile,
                                  std::string key,
                                  std::string type="string" );

        //  added @2014-11-26, to handel empty value cases, while for boolean values, empty means false
        bool Has_Value( std::string infile, std::string key, std::string type );
        bool Is_Empty( std::string line );

        //  convert string to int / double
        int     String_to_Int( std::string str );
        double  String_to_Double( std::string str );

        //    int to string
        std::string Int_to_String( int a );
        std::string IntToString( int a );
        std::string DoubleToString( double par, int precision=5 );  //  convert doubles to strings.

        void ToUpperCase( std::string &str ); //  convert to upper case
        void ToLowerCase( std::string &str ); //  convert to lower case

        bool SameStrings( std::string s1, std::string s2 );   //  compare two strings are same not not

        //  read in one string/int/double from given line by reference
        //  NOTE: these functions should be called when there is only ONE value on the right-side of '='
        void Read_String_from_Line( std::string line, std::string &value );
        void Read_Bool_from_Lile( std::string line, bool &value );
        void Read_Int_from_Line( std::string line, int &value );
        void Read_Double_from_Line( std::string line, double &value );

        //  read in one string/int/double from returned value, also, must be called when there is only ONE value
        std::string Read_String_from_Line( std::string line );
        int         Read_Int_from_Line( std::string line );
        double      Read_Double_from_Line( std::string line );

        //  for array case, I recommend to use the following
        std::string*    Read_Array_of_String_from_Line( std::string line, int &size );
        int*            Read_Array_of_Int_from_Line( std::string line, int &size );
        double*         Read_Array_of_Double_from_Line( std::string line, int &size );

        //  also, I keep these versions, but for [[known size]] of array
        void Read_Array_of_String_from_Line( std::string line,
                                             std::string str_array[],
                                             int         size,
                                             bool        warn=false );

        void Read_Array_of_Int_from_Line( std::string  line,
                                          int          int_array[],
                                          int          size,
                                          bool         warn=false );

        void Read_Array_of_Double_from_Line( std::string  line,
                                             double       double_array[],
                                             int          size,
                                             bool         warn=false );

        //  ============================================================================================
        //  read int a single value
        void Read_Value_from_File( std::string infile, std::string key, std::string& value );
        void Read_Value_from_File( std::string infile, std::string key, int& value );
        void Read_Value_from_File( std::string infile, std::string key, double& value );

        //  BRead_*_*_*
        //  'B' means that these functions will return a bool value
        bool BRead_Value_from_File( std::string infile, std::string key, std::string &value );
        bool BRead_Value_from_File( std::string infile, std::string key, int &value );
        bool BRead_Value_from_File( std::string infile, std::string key, double &value );

        std::string Read_String_from_File(std::string infile, std::string key );
        bool        Read_Bool_from_File( std::string infile, std::string key );
        int         Read_Int_from_File( std::string infile, std::string key );
        double      Read_Double_from_File( std::string infile, std::string key );

        //  read in an array of values, and the number of readed values will be also returned
        std::string* Read_Array_of_String_from_File( std::string infile,
                                                     std::string key,
                                                     int         &array_size );

        int*         Read_Array_of_Int_from_File( std::string infile,
                                                  std::string key,
                                                  int         &array_size );

        double*      Read_Array_of_Double_from_File( std::string infile,
                                                     std::string key,
                                                     int         &array_size );

        //  for know size
        void Read_Array_from_File( std::string infile,
                                   std::string key,
                                   std::string array[],
                                   int         array_size );

        void Read_Array_from_File( std::string infile,
                                   std::string key,
                                   int         array[],
                                   int         array_size);

        void Read_Array_from_File( std::string infile,
                                   std::string key,
                                   double      array[],
                                   int         array_size );

        //  modified version for given size of array
        void Read_Array_of_String_from_File(  std::string infile,
                                              std::string key,
                                              std::string str_array[],
                                              int         array_size );

        void Read_Array_of_Int_from_File( std::string infile,
                                          std::string key,
                                          int         int_array[],
                                          int         array_size );

        void Read_Array_of_Double_from_File( std::string infile,
                                             std::string key,
                                             double      double_array[],
                                             int         array_size );

        //    BRead_*_*_*
        //    B means that these functions will return a bool value  (false means failed to read array)
        bool BRead_Array_from_File( std::string infile,
                                    std::string key,
                                    std::string array[],
                                    int         array_size );

        bool BRead_Array_from_File( std::string infile,
                                    std::string key,
                                    int         array[],
                                    int         array_size );

        bool BRead_Array_from_File( std::string infile,
                                    std::string key,
                                    double      array[],
                                    int         array_size );
    }

    namespace Info{
        void WarningInfo( std::string warninginfo, bool stop=false );
        void ErrorInfo( std::string errorinfo );
    }
}

}

#endif    // __IMCMC_PARSER__
