//  TODO: the new parser class IniFile has been declared in header 'parser++.hpp', but it has not been implemented yet

#include "parser++.hpp"

namespace imcmc{
namespace parser{
    namespace Read{

    //  ======================================================
    //  test whether a file exist, added @2014-11-27
        bool Has_File( std::string paramfile ){
            bool has = false;
            std::ifstream infile(paramfile.c_str());

            if( infile ){
                has = true;
                infile.close();
            }

            return has;
        }

    //  =======================================================
    //  convert to upper case
        void ToUpperCase( std::string &str ){
            for(std::string::size_type i=0; i!=str.size(); ++i)
                str[i] = toupper(str[i]);
        }

    //  =======================================================
    //  convert to lower case
        void ToLowerCase( std::string &str ){
            for(std::string::size_type i=0; i!=str.size(); ++i)
                str[i] = tolower(str[i]);
        }

    //  =======================================================
    //  compare two strings, the case is irrelavent
        bool SameStrings( std::string s1, std::string s2 ){
            bool same = false;
            std::string S1(s1);
            std::string S2(s2);
            ToUpperCase(S1);
            ToUpperCase(S2);

            if( S1 == S2 )
                same = true;

            return same;
        }

    //  =======================================================
    //  Determine whether the line is commented out
        bool Is_Commented( std::string line ){
            bool is_commented = false;
            const std::string comments("#;");
            std::string::size_type idx = line.find_first_not_of(" \t");

        //    @2015-4-30: found a bug here, I should first check that idx is NOT std::string::npos, which
        //    amounts to check whether line is empty or not.
            if( idx != std::string::npos ){
                //  first non-empty character should not one of "#;"
                if( comments.find(line[idx]) != std::string::npos )
                    is_commented = true;
            }

            return is_commented;
        }

    //  =====================================================================
    //  @2014-11-26
    //  Need to add a new function to extract un-commented parts of a line
    //  Note: this function should be applied only to NON-commented lines
        std::string Remove_Comments( std::string line ){

           if( Is_Commented(line) ){
               // Info::ErrorInfo( "the line: " + line + " is commented out !");
             MPI_IMCMC_ERROR( "the line: " + line + " is commented out !");
           }
            
            // StopOnError( Is_Commented(line), "the line: " + line + " is commented out !" );

            std::string::size_type idx_of_comments = line.find_first_of("#;");

            return line.substr(0, idx_of_comments);
        }

        void Remove_TabSpace( std::string& value, bool info ){
            int n_space = 0;
            int n_tab = 0;

            std::string new_value = "";

            for( std::string::size_type idx = 0; idx != value.size(); ++idx ){
                if( value[idx] == ' ' )
                    ++n_space;
                else if( value[idx] == '\t' )
                    ++n_tab;
                else
                    new_value.push_back(value[idx]);
            }

            if( info ){
                std::cout << "--- imcmc::parser::Remove_TabSpace() ==>\n"
                        << "--- removed " << n_space << " space\n"
                        << "--- removed " << n_tab << " tabs\n";
            }

            value = new_value;
        }

        bool Is_Empty( std::string line ){
            bool empty = false;

            if( line.empty() ){
            //  case 1) line is truly an empty string, has no emelent(s)
                empty = true;
            }
            else if( line.find_first_not_of(" ") == std::string::npos ){
            //  case 2) line is not empty, but all elements are space(s) or blank
                empty = true;
            }

            return empty;
        }

    //  =======================================================================
    //  A(:) means that A is an array, but its size is unknown ...
    //  this is not used for the moment ...
        bool Is_Array( std::string param ){
            bool is_array = false;
            if( param.find("(:)") != std::string::npos )
                is_array = true;

            return is_array;
        }

        std::string GetArrayName( std::string param ){
            std::string::size_type idx = param.find("(:)");
            return param.substr(0, idx);
        }

    //  ========================================================================
    //  get the begin position of value
        std::string::size_type Begin_of_Value( std::string line ){

            std::string::size_type idx = line.find_first_of("=") + 1;

            if( idx == std::string::npos ){
                // Info::WarningInfo( "no \'=\' found" );
                MPI_IMCMC_WARNING( "no \'=\' found" );
                return std::string::npos;
            }
            else{   //  values are sperated by spaces or comma
                while( line[idx] == ' ' || line[idx] == ',' || line[idx] == '\t' ) ++idx;
            }

            return idx;
        }

    //  ######################################################################
        std::string::size_type End_of_Value( std::string line ){
            std::string::size_type idx_begin_of_value   = Begin_of_Value(line);
            std::string::size_type idx_end_of_value     = line.size()-1;
            std::string::size_type idx_of_comments      = line.find_first_of("#;");

            if( idx_of_comments != std::string::npos ){
                idx_end_of_value = idx_of_comments - 1;
            }
            else{
                while(  line[idx_end_of_value] == ' ' ||
                        line[idx_end_of_value] == ',' ||
                        line[idx_end_of_value] == '\t') {
                    --idx_end_of_value;
                }
            }

            if( idx_end_of_value < idx_begin_of_value )
                return std::string::npos;
            else
                return idx_end_of_value;
        }
    //  ######################################################################

    //  Check whether the key is included in line, if correct, one key should be found.
        bool Has_Key( std::string line, std::string key ){
            bool has_key = true;

            std::string::size_type idx = line.find_first_of("=");
            std::string subline = line.substr(0, idx);
            std::istringstream stream(subline);
            std::string temp_key;

            std::string copy_of_key = key;
            ToUpperCase(copy_of_key);

            int nkey=0;

            while(stream >> temp_key){
                ToUpperCase(temp_key);
                if( temp_key == copy_of_key )
                    ++nkey;
            }
            
            // StopOnError(nkey>1,"more than one " + key + " ...");

            if(nkey < 1)
                has_key = false;
            else if(nkey == 1)
                has_key = true;
           else if(nkey > 1){
               // Info::ErrorInfo( "more than one " + key + " ..." );
                MPI_IMCMC_ERROR( "more than one " + key + " ..." );
           }

            return has_key;
        }

        bool Has_Key_in_File( std::string file, std::string key ){

            bool has = false;
            std::ifstream infile(file.c_str());
            std::string line;
            int count = 0;

            while( std::getline(infile, line) ){
                if( !Is_Commented( line ) && Has_Key(line, key) )
                    ++count;
            }

        //    @2015-4-30, a bug was found here
            infile.close();    //    I forgot to add this line ...

            // StopOnError( count>1, "more than one " + key + " in paramfile, ambiguous" );

            if( count == 0 ){
                has = false;
            }
            else if( count == 1 ){
                has = true;
            }
           else{
               // Info::ErrorInfo( "more than one " + key + " in paramfile, ambiguous" );
                MPI_IMCMC_ERROR( "more than one " + key + " in paramfile, ambiguous" );
           }
            return has;
        }

    //  number of values
        int Num_of_Value(   std::string line,
                            std::string::size_type begin,
                            std::string::size_type end,
                            std::string type ){

            int num_of_value    = 0;       //  at least one value
            std::string subline = line.substr( begin, end-begin+1 );

            for( std::string::size_type i=0; i!=subline.size(); ++i ){
                if( subline[i] == ',' || subline[i] == '\t')    //    @2015-4-30
                    subline[i] = ' ';
            }

            std::istringstream stream(subline);

            if( SameStrings(type,"string") ){
                std::string value;
                while( stream >> value ) ++num_of_value;
            }
            else if( SameStrings(type,"int") ){
                int value;
                while( stream >> value ) ++num_of_value;
            }
            else if( SameStrings(type,"double") ){
                double value;
                while( stream >> value ) ++num_of_value;
            }
            // else{
            //     StopOnError( true, "unrecognized data type" );
            // }
           else{
               // Info::ErrorInfo( "unrecognized data type" );
                MPI_IMCMC_ERROR( "unrecognized data type" );
           }

            return num_of_value;
        }

        int Num_of_Value_for_Key(   std::string infile,
                                    std::string key,
                                    std::string type    ){

            std::string line = Read_String_from_File(infile, key);
            std::string::size_type idx_begin = Begin_of_Value(line);
            std::string::size_type idx_end = End_of_Value(line);

            return Num_of_Value(line, idx_begin, idx_end, type );
        }

    //  added @2014-11-26
        bool Has_Value( std::string infile, std::string key, std::string type ){

            int value_num = 0;

            if( Has_Key_in_File( infile, key ) )
                value_num = Num_of_Value_for_Key( infile, key, type );

            if( value_num > 0 )
                return true;
            else
                return false;
        }

        int String_to_Int( std::string str ){
            std::istringstream stream(str);
            int value;
            stream >> value;    //    this makes sure that only the first value is extracted
            return value;
//            return std::stoi(str);    //  C++11, add -std=c++11 to CFLAG
        }

        double String_to_Double( std::string str ){
            std::istringstream stream(str);
            double value;
            stream >> value;    //    this makes sure that only the first value is extracted
            return value;
//            return std::stod(str);  //  C++11, add -std=c++11 to CFLAG
        }

        std::string DoubleToString( double par, int precision ){    //    default precision is 5 digits after decimal
            std::stringstream stream;
            stream << std::setprecision(precision) << par;

            std::string par_in_string;
            stream >> par_in_string;

            stream.clear();

            return par_in_string;
        }

    //    added @2014-11-28
        std::string Int_to_String( int a ){
            std::stringstream stream;
            stream << a;
            return stream.str();
        }

        std::string IntToString( int a ){
            std::stringstream stream;
            stream << a;
            return stream.str();
        }

    //  --------------------------------------------------------------------------------------------
        void Read_String_from_Line( std::string line, std::string &value ){
        //  bug-fix: Mar-9-2015
        //  before finding the begin & end of values, remove the inline comments first.
        //  Remove_Comments() has been there for a very long time, but I never used it, what a shame...
            std::string line_temp = Remove_Comments(line);
            std::string::size_type idx_begin = Begin_of_Value(line_temp);
            std::string::size_type idx_end = End_of_Value(line_temp);
            value = line_temp.substr(idx_begin, idx_end-idx_begin+1);
        }

        void Read_Int_from_Line( std::string line, int &value ){
            std::string temp;
            Read_String_from_Line( line, temp );
            value = String_to_Int(temp);
        }

        void Read_Double_from_Line( std::string line, double &value ){
            std::string temp;
            Read_String_from_Line( line, temp );
            value = String_to_Double(temp);
        }

    //  other versions ...
        std::string Read_String_from_Line( std::string line ){
//            std::string::size_type idx_begin = Begin_of_Value(line);
//            std::string::size_type idx_end = End_of_Value(line);

        //  bug-fix: Mar-9-2015
            std::string line_temp = Remove_Comments(line);
            std::string::size_type idx_begin = Begin_of_Value(line_temp);
            std::string::size_type idx_end = End_of_Value(line_temp);
            return line.substr(idx_begin, idx_end-idx_begin+1);
        }

        int Read_Int_from_Line( std::string line ){
            std::string temp;
            Read_String_from_Line( line, temp );
            return String_to_Int(temp);
        }

        double Read_Double_from_Line( std::string line ){
            std::string temp;
            Read_String_from_Line( line, temp );
            return String_to_Double(temp);
        }

//        int Get_Array_of_String( std::string line, std::string *str_array ){
        std::string* Read_Array_of_String_from_Line( std::string line, int &size ){
            int num_of_value = 0;
            std::vector<std::string> array;
            std::string::size_type idx_begin = Begin_of_Value(line);
            std::string::size_type idx_end = End_of_Value(line);
            std::string subline = line.substr( idx_begin, idx_end-idx_begin+1 );

            for( std::string::size_type i=0; i!=subline.size(); ++i ){
                if( subline[i] == ',' ) subline[i] = ' ';
            }

            std::istringstream stream(subline);
            std::string value;

            while( stream >> value ){
                array.push_back(value);
                ++num_of_value;
            }

            size = num_of_value;

            std::string *str_array = new std::string[num_of_value];  // remember to delete[]
            for(int i=0; i<num_of_value; ++i) str_array[i] = array[i];

            return str_array;
        }

        int* Read_Array_of_Int_from_Line( std::string line, int &size  ){

            int num_of_value = 0;
            std::vector<int> array;

            std::string::size_type idx_begin = Begin_of_Value(line);
            std::string::size_type idx_end = End_of_Value(line);
            std::string subline = line.substr( idx_begin, idx_end-idx_begin+1 );

            for( std::string::size_type i=0; i!=subline.size(); ++i ){
                if( subline[i] == ',' ) subline[i] = ' ';
            }

            std::istringstream stream(subline);
            int value;

            while( stream >> value ){
                array.push_back(value);
                ++num_of_value;
            }

            size = num_of_value;
            int *value_array = new int[num_of_value];  // remember to delete[]
            for(int i=0; i<num_of_value; ++i) value_array[i] = array[i];

            return value_array;
        }

        double *Read_Array_of_Double_from_Line( std::string line, int &size ){

            int num_of_value = 0;
            std::vector<double> array;

            std::string::size_type idx_begin = Begin_of_Value(line);
            std::string::size_type idx_end = End_of_Value(line);
            std::string subline = line.substr( idx_begin, idx_end-idx_begin+1 );

            for( std::string::size_type i=0; i!=subline.size(); ++i ){
                if( subline[i] == ',' ) subline[i] = ' ';
            }

            std::istringstream stream(subline);
            double value;

            while( stream >> value ){
                array.push_back(value);
                ++num_of_value;
            }

            size = num_of_value;

            double *value_array = new double[num_of_value];  // remember to delete[]
            for(int i=0; i<num_of_value; ++i) value_array[i] = array[i];

            return value_array;
        }

        void Read_Array_of_String_from_Line(    std::string line,
                                                std::string str_array[],
                                                int size,
                                                bool warn   ){

            int num_of_value = 0;
            std::vector<std::string> array;
            std::string::size_type idx_begin = Begin_of_Value(line);
            std::string::size_type idx_end = End_of_Value(line);
            std::string subline = line.substr( idx_begin, idx_end-idx_begin+1 );

            for( std::string::size_type i=0; i!=subline.size(); ++i ){
                if( subline[i] == ',' ) subline[i] = ' ';
            }

            std::istringstream stream(subline);
            std::string value;

            while( stream >> value ){
                array.push_back(value);
                ++num_of_value;
            }

            if( size == num_of_value ){
                for(int i=0; i<num_of_value; ++i)
                    str_array[i] = array[i];
            }
            else if(warn){
                // Info::WarningInfo( "number of values does not match" );
                MPI_IMCMC_WARNING( "number of values does not match" );
            }
        }

        void Read_Array_of_Int_from_Line(   std::string line,
                                            int int_array[],
                                            int size,
                                            bool warn   ){

            int num_of_value = 0;
            std::vector<int> array;

            std::string::size_type idx_begin = Begin_of_Value(line);
            std::string::size_type idx_end = End_of_Value(line);
            std::string subline = line.substr( idx_begin, idx_end-idx_begin+1 );

            for( std::string::size_type i=0; i!=subline.size(); ++i ){
                if( subline[i] == ',' ) subline[i] = ' ';
            }

            std::istringstream stream(subline);
            int value;

            while( stream >> value ){
                array.push_back(value);
                ++num_of_value;
            }

            if( size == num_of_value ){
                for(int i=0; i<num_of_value; ++i)
                    int_array[i] = array[i];
            }
            else if(warn){
                // Info::WarningInfo( "number of values does not match" );
                MPI_IMCMC_WARNING( "number of values does not match" );
            }
        }

        void Read_Array_of_Double_from_Line(    std::string line,
                                                double double_array[],
                                                int size,
                                                bool warn   ){

            int num_of_value = 0;
            std::vector<double> array;

            std::string::size_type idx_begin = Begin_of_Value(line);
            std::string::size_type idx_end = End_of_Value(line);
            std::string subline = line.substr( idx_begin, idx_end-idx_begin+1 );

            for( std::string::size_type i=0; i!=subline.size(); ++i ){
                if( subline[i] == ',' ) subline[i] = ' ';
            }

            std::istringstream stream(subline);
            double value;

            while( stream >> value ){
                array.push_back(value);
                ++num_of_value;
            }

            if( size == num_of_value ){
                for(int i=0; i<num_of_value; ++i)
                    double_array[i] = array[i];
            }
            else if(warn){
                // Info::WarningInfo( "number of values does not match" );
                MPI_IMCMC_WARNING( "number of values does not match" );
            }
        }

    //  ====================================
    //      Key-Value ( single value )
    //  ====================================

    //  Value will be stored as string, this is the default case.
    //  Luckily, thanks to the overloading ability of C++, we can use the same function name to do
    //  many other things
        void Read_Value_from_File(  std::string infile,
                                    std::string key,
                                    std::string &value  ){

            bool readed = false;
            std::string line;
            std::ifstream in(infile.c_str());

            while( std::getline(in, line) ){
                if( !Is_Commented( line ) && Has_Key(line, key) && Has_Value(infile, key, "string") ){   // found key
                    Read_String_from_Line(line, value);
                    readed = true;
                }
            }

            in.close();
            
            // StopOnError( readed==false, "Failed to read value for " + key );
           if( !readed ){
               // Info::ErrorInfo( "Failed to read value for " + key );
                MPI_IMCMC_ERROR( "Failed to read value for " + key );
           }
        }

    //  Now some overloaded functions of the above Read_Value
        void Read_Value_from_File(  std::string infile,
                                    std::string key,
                                    int &value  ){

            bool readed = false;
            std::string line;
            std::ifstream in(infile.c_str());

            while( std::getline(in, line) ){
                if( !Is_Commented( line ) && Has_Key(line, key) && Has_Value(infile, key, "int") ){   // found key
                    Read_Int_from_Line(line, value);
                    readed = true;
                }
            }

            in.close();
           if( !readed ){
               // Info::ErrorInfo( "Failed to read value for " + key );
                MPI_IMCMC_ERROR( "Failed to read value for " + key );
           }
            // StopOnError( readed==false, "Failed to read value for " + key );
        }

        void Read_Value_from_File(  std::string infile,
                                    std::string key,
                                    double &value   ){

            bool readed = false;
            std::string line;
            std::ifstream in(infile.c_str());

            while( std::getline(in, line) ){
                if( !Is_Commented( line ) && Has_Key(line, key) && Has_Value(infile, key, "double") ){   // found key
                    Read_Double_from_Line(line, value);
                    readed = true;
                }
            }

            in.close();
           if( !readed ){
               // Info::ErrorInfo( "Failed to read value for " + key );
                MPI_IMCMC_ERROR( "Failed to read value for " + key );
           }
        }

    //  bool version
        bool BRead_Value_from_File( std::string infile,
                                    std::string key,
                                    std::string &value ){

            bool readed = false;
            std::string line;
            std::ifstream in(infile.c_str());

            while( std::getline(in, line) ){
                if( !Is_Commented( line ) && Has_Key(line, key) && Has_Value(infile, key, "string") ){   // found key
                    Read_String_from_Line(line, value);
                    readed = true;
                }
            }

            in.close();
            return readed;
        }

    //  Now some overloaded functions of the above Read_Value
        bool BRead_Value_from_File( std::string infile, std::string key, int &value ){

            bool readed = false;
            std::string line;
            std::ifstream in(infile.c_str());

            while( std::getline(in, line) ){
                if( !Is_Commented( line ) && Has_Key(line, key) && Has_Value(infile, key, "int") ){   // found key
                    Read_Int_from_Line(line, value);
                    readed = true;
                }
            }

            in.close();
            return readed;
        }

        bool BRead_Value_from_File( std::string infile, std::string key, double &value ){

            bool readed = false;
            std::string line;
            std::ifstream in(infile.c_str());

            while( std::getline(in, line) ){
                if( !Is_Commented( line ) && Has_Key(line, key) && Has_Value(infile, key, "double") ){   // found key
                    Read_Double_from_Line(line, value);
                    readed = true;
                }
            }

            in.close();
            return readed;
        }

        std::string Read_String_from_File( std::string infile, std::string key ){

            bool readed = false;
            std::string value;
            std::string line;
            std::ifstream in(infile.c_str());

            while( std::getline(in, line) ){
                if( !Is_Commented( line ) && Has_Key(line, key) ){   // found key
                    Read_String_from_Line(line, value);
                    readed = true;
                }
            }

            in.close();
           if( !readed ){
               // Info::ErrorInfo( "Failed to read value for " + key );
                MPI_IMCMC_ERROR( "Failed to read value for " + key );
           }
            // StopOnError( readed==false, "Failed to read value for " + key );
            
            return value;
        }

        bool Read_Bool_from_File( std::string infile, std::string key ){

            bool readed = false;
            bool bool_value = false;
            std::string value;
            std::string line;
            std::ifstream in(infile.c_str());

            if( Has_Key_in_File(infile, key) ){
                while( std::getline(in, line) ){
                    if( !Is_Commented( line ) && Has_Key(line, key) ){   // found key
                        Read_String_from_Line(line, value);
                        readed = true;
                    }
                }

                in.close();

                if( !readed ){
                    // Info::WarningInfo( "Failed to read value for "+key+", so false is return");
                    MPI_IMCMC_WARNING( "Failed to read value for "+key+", so false is return");
                    bool_value = false;
                }
                else{
                    ToUpperCase(value);

                    Remove_TabSpace(value);

                    if( value == "T" || value == "TRUE" )
                        bool_value = true;
                    else if( value == "F" || value == "FALSE" || value == "" )
                        bool_value = false;
                    else{
                       // Info::ErrorInfo( "unrecognized bool value: " + value );
                        MPI_IMCMC_ERROR( "unrecognized bool value: " + value );
                        // StopOnError( true, "unrecognized bool value: " + value );
                    }
                }
            }

            return bool_value;
        }

        int Read_Int_from_File( std::string infile, std::string key ){

            bool readed = false;
            int value;
            std::string line;
            std::ifstream in(infile.c_str());

            while( std::getline(in, line) ){
                if( !Is_Commented( line ) && Has_Key(line, key) && Has_Value(infile, key, "int") ){   // found key
                    Read_Int_from_Line(line, value);
                    readed = true;
                }
            }

            in.close();
           if( !readed ){
               // Info::ErrorInfo( "Failed to read value for " + key );
                MPI_IMCMC_ERROR( "Failed to read value for " + key );
           }
            // StopOnError( readed==false, "Failed to read value for " + key );
            return value;
        }

        double Read_Double_from_File( std::string infile, std::string key ){

            bool readed = false;
            double value;
            std::string line;
            std::ifstream in(infile.c_str());

            while( std::getline(in, line) ){
                if( !Is_Commented( line ) && Has_Key(line, key) && Has_Value(infile, key, "double") ){   // found key
                    Read_Double_from_Line(line, value);
                    readed = true;
                }
            }

            in.close();
           if( !readed ){
               // Info::ErrorInfo( "Failed to read value for " + key );
                MPI_IMCMC_ERROR( "Failed to read value for " + key );
           }
            // StopOnError( readed==false, "Failed to read value for " + key );
            return value;
        }

    //  ===================================
    //  Key-Value ( array of values )
    //  ===================================
        std::string* Read_Array_of_String_from_File( std::string infile, std::string key, int &array_size ){

            bool readed = false;
            std::string line;
            std::string *array = NULL;
            std::ifstream in(infile.c_str());

            while( std::getline(in, line) ){
                if( !Is_Commented( line ) && Has_Key(line, key) ){   // found key
                    array = Read_Array_of_String_from_Line( line, array_size );
                    readed = true;
                }
            }

            in.close();
            if( !readed ){
               // Info::WarningInfo("Failed to read value for "+key);
                MPI_IMCMC_WARNING("Failed to read value for "+key);
                // DetectWarning( readed==false, "Failed to read value for "+key );
                return NULL;
            }
            else
                return array;
        }

        int* Read_Array_of_Int_from_File( std::string infile, std::string key, int &array_size ){

            bool readed = false;
            std::string line;
            int *array = NULL;
            std::ifstream in(infile.c_str());

            while( std::getline(in, line) ){
                if( !Is_Commented( line ) && Has_Key(line, key) ){   // found key
                    array = Read_Array_of_Int_from_Line( line, array_size );
                    readed = true;
                }
            }

            in.close();
            if( !readed ){
                // Info::WarningInfo("Failed to read value for "+key);
                MPI_IMCMC_WARNING("Failed to read value for "+key);
                return NULL;
            }
            else
                return array;
        }

        double* Read_Array_of_Double_from_File( std::string infile, std::string key, int &array_size ){

            bool readed = false;
            std::string line;
            double *array = NULL;
            std::ifstream in(infile.c_str());

            while( std::getline(in, line) ){
                if( !Is_Commented( line ) && Has_Key(line, key) ){   // found key
                    array = Read_Array_of_Double_from_Line( line, array_size );
                    readed = true;
                }
            }

            in.close();
            if( !readed ){
                // Info::WarningInfo("Failed to read value for "+key);
                MPI_IMCMC_WARNING("Failed to read value for "+key);
                return NULL;
            }
            else
                return array;
        }

    //  the following are different version, for known size of array
        void Read_Array_of_String_from_File( std::string infile, std::string key, std::string str_array[], int array_size ){

            bool readed = false;
            std::string line;
            std::ifstream in(infile.c_str());

            while( std::getline(in, line) ){
                if( !Is_Commented( line ) && Has_Key(line, key) ){   // found key
                    Read_Array_of_String_from_Line( line, str_array, array_size );
                    readed = true;
                }
            }
            in.close();
           if( !readed ){
               // Info::ErrorInfo( "Failed to read str_array[] for " + key );
                MPI_IMCMC_ERROR( "Failed to read str_array[] for " + key );
           }
            // StopOnError( readed==false, "Failed to read value for " + key );
        }

        void Read_Array_of_Int_from_File( std::string infile, std::string key, int int_array[], int array_size ){

            bool readed = false;
            std::string line;
            std::ifstream in(infile.c_str());

            while( std::getline(in, line) ){
                if( !Is_Commented( line ) && Has_Key(line, key) ){   // found key
                    Read_Array_of_Int_from_Line( line, int_array, array_size );
                    readed = true;
                }
            }

            in.close();
           if( !readed ){
               // Info::ErrorInfo( "Failed to read int_array[] for " + key );
                MPI_IMCMC_ERROR( "Failed to read int_array[] for " + key );
           }
            // StopOnError( readed==false, "Failed to read value for " + key );
        }

        void Read_Array_of_Double_from_File( std::string infile, std::string key, double double_array[], int array_size ){

            bool readed = false;
            std::string line;
//            double *array;
            std::ifstream in(infile.c_str());

            while( std::getline(in, line) ){
                if( !Is_Commented( line ) && Has_Key(line, key) ){   // found key
                    Read_Array_of_Double_from_Line( line, double_array, array_size );
                    readed = true;
                }
            }

            in.close();
           if( !readed ){
               // Info::ErrorInfo( "Failed to read double_array[] for " + key );
                MPI_IMCMC_ERROR( "Failed to read double_array[] for " + key );
           }
            // StopOnError( readed==false, "Failed to read value for " + key );
        }

        void Read_Array_from_File(std::string infile, std::string key, std::string array[], int array_size){
            Read_Array_of_String_from_File( infile, key, array, array_size );
        }

        void Read_Array_from_File(std::string infile, std::string key, int array[], int array_size){
            Read_Array_of_Int_from_File( infile, key, array, array_size );
        }

        void Read_Array_from_File(std::string infile, std::string key, double array[], int array_size){
            Read_Array_of_Double_from_File( infile, key, array, array_size );
        }


        bool BRead_Array_from_File(std::string infile, std::string key, std::string array[], int array_size){

            bool readed = false;
            std::string line;
            std::ifstream in(infile.c_str());

            while( std::getline(in, line) ){
                if( !Is_Commented( line ) && Has_Key(line, key) ){   // found key
                    Read_Array_of_String_from_Line( line, array, array_size );
                    readed = true;
                }
            }
            in.close();
            return readed;
        }

        bool BRead_Array_from_File(std::string infile, std::string key, int array[], int array_size){

            bool readed = false;
            std::string line;
            std::ifstream in(infile.c_str());

            while( std::getline(in, line) ){
                if( !Is_Commented( line ) && Has_Key(line, key) ){   // found key
                    Read_Array_of_Int_from_Line( line, array, array_size );
                    readed = true;
                }
            }
            in.close();
            return readed;
        }

        bool BRead_Array_from_File(std::string infile, std::string key, double array[], int array_size){

            bool readed = false;
            std::string line;
            std::ifstream in(infile.c_str());

            while( std::getline(in, line) ){
                if( !Is_Commented( line ) && Has_Key(line, key) ){   // found key
                    Read_Array_of_Double_from_Line( line, array, array_size );
                    readed = true;
                }
            }
            in.close();
            return readed;
        }

    }

    namespace Info{
    //  Though imcmc has two functions similar to the follwings, I still keep them here, because these
    //  two are used inside parser only, and NO imcmc function were used here, so this parser is actually
    //  completely independent and can be used anywhere else directly (don't forget the header)

    //  print warning information, but will not stop execution unless stop=true
        void WarningInfo( std::string warninginfo, bool stop ){

            std::cout << "\n#########################################################\n"
                      << "#--- File Name: " << __FILE__ << "\n"
                      << "#--- Line    #: " << __LINE__ << "\n"
                      << "#--- Func Name: " << __FUNCTION__ << "\n"
                      << "#--- Warning Infomation: " << warninginfo << "\n\n";

            if(stop){
                std::cout << "---->\tExit from " << __FUNCTION__ << "\n";
                throw std::runtime_error( "you requested to stop here!\n" );
            }
        }

    //  print error information and stop
        void ErrorInfo( std::string errorinfo ){

            std::cout << "\n#########################################################\n"
                      << "#--- File Name: " << __FILE__ << "\n"
                      << "#--- Line    #: " << __LINE__ << "\n"
                      << "#--- Func Name: " << __FUNCTION__ << "\n";

            throw std::runtime_error( errorinfo + "\n" );
        }
    }
}
}
