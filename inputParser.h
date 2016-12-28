/*
 * input parser inspired by
 * http://stackoverflow.com/questions/865668/how-to-parse-command-line-arguments-in-c
 * 
 */

#ifndef INPUT_PARSER_H_LBDEM
#define INPUT_PARSER_H_LBDEM

#include <algorithm>
#include <vector>
#include <string>
#include <exception>
#include <stdexcept>

class InputParser{
public:
  InputParser (int &argc, char **argv)
  {
    for (int i=1; i < argc; ++i)
      this->tokens.push_back(std::string(argv[i]));
  }
  /// @author iain

  int readInt(const char* option) const
  {
    std::string strval = getCmdOption(option);
    return atoi(strval.c_str());
  }

  double readDouble(const char* option) const
  {
    std::string strval = getCmdOption(option);
    return atof(strval.c_str());
  }

  bool readBool(const char* option) const
  {
    std::string strval = getCmdOption(option);
    return (strval != "0") && (strval != "false");
  }

  std::string readString(const char* option) const
  {
    std::string strval = getCmdOption(option);
    return strval;    
  }
  const std::string getCmdOption(const std::string &option) const
  {
    std::vector<std::string>::const_iterator itr;
    itr =  std::find(this->tokens.begin(), this->tokens.end(), option);
    if (itr != this->tokens.end() && ++itr != this->tokens.end()){
      return *itr;
    }
    throw std::invalid_argument(option);
    return "";
  }
  /// @author iain
  bool cmdOptionExists(const std::string &option) const
  {
    return std::find(this->tokens.begin(), this->tokens.end(), option)
      != this->tokens.end();
  }
private:
  std::vector <std::string> tokens;
};

#endif
