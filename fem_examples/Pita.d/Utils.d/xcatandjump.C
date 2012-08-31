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
    std::cout << "Syntax: xcatandjump filename [#files]" << std::endl;
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
  int currNumFile = 0;
  
  // All files are present and open

  int integerBuffer;
  double floatBuffer;
  double timeStamp;
  std::string line;
  std::stringstream lineStream;

  // Get first line
  std::getline(partFileStream, line, '\n');
  catFileStream << line << '\n'; 
  jumpFileStream << line << '\n';

  // Get number of dofs
  std::getline(partFileStream, line, '\n');
  catFileStream << line << '\n';
  jumpFileStream << line << '\n';
  lineStream.clear();
  lineStream.str(line);
  lineStream >> integerBuffer;
  unsigned int numDofs = integerBuffer < 0 ? 0u : integerBuffer;

  // Keep values read at last time stamp to compute the jump values
  std::vector<std::vector<double> > values(numDofs);
  std::vector<std::vector<double> >::iterator it;

  // Get time stamp 
  std::getline(partFileStream, line, '\n');
  catFileStream << line << '\n';
  lineStream.clear();
  lineStream.str(line);
  lineStream >> floatBuffer;
  timeStamp = floatBuffer;
  
  // Store values
  for (it = values.begin(); it != values.end(); ++it)
  {
    std::getline(partFileStream, line, '\n');
    catFileStream << line << '\n';
    lineStream.clear();
    lineStream.str(line);
    lineStream >> floatBuffer;
    while (lineStream)
    {
      it->push_back(floatBuffer);
      lineStream >> floatBuffer;
    }
  }

  while (true) 
  {
    // Get time stamp
    std::getline(partFileStream, line, '\n');
    if (line == "")
      break;
    catFileStream << line << '\n';
    lineStream.str(line);
    lineStream >> floatBuffer; 
    timeStamp = floatBuffer;
  
    // Store values
    for (it = values.begin(); it != values.end(); ++it)
    {
      std::getline(partFileStream, line, '\n');
      catFileStream << line << '\n';
      lineStream.clear();
      lineStream.str(line);
      lineStream >> floatBuffer;
      for (std::vector<double>::iterator jt = it->begin(); jt != it->end(); ++jt)
      {
        *jt = floatBuffer;
        lineStream >> floatBuffer;
      }
    }
  }
 
  partFileStream.close();

  while (true)
  {
    ++currNumFile;
    if (currNumFile == numFiles)
      break;
    
    std::stringstream partFileNameStream;
    partFileNameStream << baseFileName << '.' << currNumFile;
    partFileName = partFileNameStream.str();
    partFileStream.clear();
    partFileStream.open(partFileName.c_str());
    if (!partFileStream)
    {
      numFiles = currNumFile;
      break;
    }
    
    // Ignore first line
    std::getline(partFileStream, line, '\n');

    // Check number of dofs
    std::getline(partFileStream, line, '\n');
    lineStream.clear();
    lineStream.str(line);
    lineStream >> integerBuffer;
    if (numDofs != static_cast<unsigned int>(integerBuffer))
    {
      std::cerr << "Inconsistent number of dofs\n";
      return 1;
    }
    
    // Get time stamp 
    std::getline(partFileStream, line, '\n');
    catFileStream << line << '\n';
    lineStream.clear();
    lineStream.str(line);
    lineStream >> floatBuffer;
    timeStamp = floatBuffer;
    
    // Compute jump
    jumpFileStream << line << '\n';
    for (it = values.begin(); it != values.end(); ++it)
    {
      std::getline(partFileStream, line, '\n');
      catFileStream << line << '\n';
      lineStream.clear();
      lineStream.str(line);
      lineStream >> floatBuffer;
      for (std::vector<double>::iterator jt = it->begin(); jt != it->end(); ++jt)
      {
        jumpFileStream << ' ' << std::scientific << floatBuffer - *jt;
        *jt = floatBuffer;
        lineStream >> floatBuffer;
      }
      jumpFileStream << '\n';
    }

    while (true) 
    {
      // Get time stamp
      std::getline(partFileStream, line, '\n');
      if (line == "")
        break;
      catFileStream << line << '\n';
      lineStream.clear();
      lineStream.str(line);
      lineStream >> floatBuffer; 
      timeStamp = floatBuffer;
    
      // Store values
      for (it = values.begin(); it != values.end(); ++it)
      {
        std::getline(partFileStream, line, '\n');
        catFileStream << line << '\n';
        lineStream.clear();
        lineStream.str(line);
        lineStream >> floatBuffer;
        for (std::vector<double>::iterator jt = it->begin(); jt != it->end(); ++jt)
        {
          *jt = floatBuffer;
          lineStream >> floatBuffer;
        }
      }
    }

    partFileStream.clear();
    partFileStream.close();
  }
    
  jumpFileStream.close();
  catFileStream.close();

  std::cout << numFiles << " files read" << std::endl;
  return 0;
}
