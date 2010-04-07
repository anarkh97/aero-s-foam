#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <algorithm>

int main (int argc, char * argv[])
{

  if (argc < 2)
  {
    std::cout << "Syntax: catjumpclean filename [#files]" << std::endl;
    return 0;
  }

  std::string baseFileName(argv[1]);
  int numFiles = (argc < 3) ? -1 : std::max(std::atoi(argv[2]), 0);

  if (numFiles == 0)
  {
    return 0;
  }

  std::string partFileName = baseFileName;
  partFileName.append(".0");
  std::ifstream partFileStream;
  partFileStream.open(partFileName.c_str()); 
  if (!partFileStream)
  {
    std::cout << "Error: Cannot open " << partFileName << std::endl;
    return 0;
  }

  std::string catFileName = baseFileName;
  catFileName.append(".cat");
  std::ofstream catFileStream(catFileName.c_str(), std::ios::trunc);

  if (!catFileStream)
  {
    std::cerr << "Error: Cannot open " << catFileName << std::endl;
    return 0;
  }

  std::string jumpFileName = baseFileName;
  jumpFileName.append(".jumps");
  std::ofstream jumpFileStream(jumpFileName.c_str(), std::ios::trunc);

  if (!jumpFileStream)
  {
    std::cerr << "Error: Cannot open " << jumpFileName << std::endl;
    return 0;
  }

  double numericBuffer;
  std::string word;

  while (partFileStream.peek() == '#')
  {
    std::getline(partFileStream, word, '\n');
  }
  std::stringstream lineStream;
  partFileStream.get(*lineStream.rdbuf());
  lineStream >> numericBuffer; 
  jumpFileStream << std::scientific << numericBuffer;
  int numValues = 0;
  while (!lineStream.eof()) 
  {
    ++numValues;
    lineStream >> numericBuffer;
    jumpFileStream << ' ' << std::scientific << double(0.0);
  }

  std::vector<double> valueList(numValues, 0.0);
  std::vector<double>::iterator it;
  std::stringstream s;
  std::string numFileName;
 
  int currNumFile = 0;
  while (true)
  {
   
    partFileStream.seekg(0, std::ios_base::beg);
    partFileStream.clear();
    
    while (!partFileStream.eof())
    {
      if (partFileStream.peek() == '#')
      {
        std::getline(partFileStream, word, '\n');
      }
      else
      {
        std::streampos startingPos = partFileStream.tellg();
        partFileStream.get(*catFileStream.rdbuf());
        while (partFileStream.peek() == '\n')
        {
          partFileStream.ignore();
        }
        if (partFileStream.eof())
        {
          partFileStream.clear();
          partFileStream.seekg(startingPos);
          partFileStream >> numericBuffer >> numericBuffer;
          it = valueList.begin();
          while (!partFileStream.eof() && it != valueList.end())
          {
            *(it++) = numericBuffer;
            partFileStream >> numericBuffer;
          }
          jumpFileStream << std::endl;
        }
        catFileStream << std::endl;
      }
    }
    partFileStream.clear();
    partFileStream.close();

    std::remove(partFileName.c_str());

    ++currNumFile;
    if (currNumFile == numFiles) break;

    s.clear();
    s.str("");
    s << '.' << currNumFile;
    s >> numFileName;
    partFileName = baseFileName + numFileName;
    partFileStream.open(partFileName.c_str());
    if (!partFileStream)
    {
      numFiles = currNumFile;
      break;
    }

    while (partFileStream.peek() == '#')
    {
      std::getline(partFileStream, word);
    }
  
    lineStream.str("");
    lineStream.seekg(0, std::ios_base::beg);
    lineStream.seekp(0, std::ios_base::end);
    lineStream.clear();
    partFileStream.get(*lineStream.rdbuf());
    lineStream >> numericBuffer;
    jumpFileStream << std::scientific << numericBuffer;
    it = valueList.begin();
    while (!lineStream.eof() && it != valueList.end())
    {
      lineStream >> numericBuffer;
      jumpFileStream << ' ' << std::scientific << std::fabs(numericBuffer - *(it++));
    }
 
  }
  
  std::cout << numFiles << " files read" << std::endl;

  jumpFileStream.close();
  catFileStream.close();

  return 0;
}
