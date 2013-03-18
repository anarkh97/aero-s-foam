#include "FileNameInfo.h"

#include <Driver.d/GeoSource.h>

#include <sstream>

extern GeoSource *geoSource;

namespace Rom {

inline
std::string
FileNameInfo::getPrefix(const std::string &str) {
  return str.substr(0, str.find('.'));
}

FileNameInfo::FileNameInfo() :
  prefix_(getPrefix(geoSource->getCheckFileInfo()->checkfile))
{}

FileNameInfo::FileNameInfo(const std::string &prefix) :
  prefix_(prefix)
{}

std::string
FileNameInfo::basisFileName(const BasisId &id) const {
  std::ostringstream builder;

 if(domain->solInfo().svdPodRom) {
  if (id.level() == 0) {
     //filePrint(stderr,"*** Performing SVD decomposition on %s \n", domain->solInfo().snapfiPodRom);
     builder << domain->solInfo().snapfiPodRom; }
  else if (id.level() == 1){
     //filePrint(stderr,"*** Saving Basis to file %s \n", domain->solInfo().SVDoutput);
     builder << domain->solInfo().SVDoutput; } }
 else {
  if(id.level() == 0) {
    if(id.type() == 0) 
        builder << domain->solInfo().statePodRomFile;
    if(id.type() == 1) 
        builder << domain->solInfo().residualPodRomFile;
    if(id.type() == 2)
        builder << domain->solInfo().jacobianPodRomFile;
    if(id.type() == 3)
        builder << domain->solInfo().forcePodRomFile;
    if(id.type() == 4)
        builder << domain->solInfo().accelPodRomFile;
    if(id.type() == 5)
        builder << domain->solInfo().velocPodRomFile;
    if(id.type() == 6)
        builder << domain->solInfo().isvPodRomFile;}
  else if(id.level() == 1) {
        builder << domain->solInfo().readInROBorModes;}
     }

  std::string mystrg;
  mystrg = builder.str();
//  filePrint(stderr,"id.level = %d, id.type = %d,  %s \n",id.level(), id.type(), mystrg.c_str());	

/*  else {
	const char suffix[] = "rob";
	builder << prefix()
          << "." << toString(id.type())
          << "." << toString(id.level())
          << "." << suffix; }*/

  return builder.str();
}

} /* end namespace Rom */
