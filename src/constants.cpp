#include "constants.hpp"

#include <json.hpp>
#include <fstream>
#include <iostream>

Constants::Constants() {}

Constants::Constants(const std::string &input) {
  readJson(input);
}

Constants::Constants(int nZ_in, int nN_in, double initialDt_in, double Ra_in, double Pr_in, int aspectRatio_in, double timeBetweenSaves_in, double totalTime_in, std::string &saveFolder_in, std::string &icFile_in)
  : nZ {nZ_in}
  , nN {nN_in}
  , initialDt {initialDt_in}
  , Ra {Ra_in}
  , Pr {Pr_in}
  , aspectRatio {aspectRatio_in}
  , timeBetweenSaves {timeBetweenSaves_in}
  , totalTime {totalTime_in}
  , saveFolder {saveFolder_in}
  , icFile {icFile_in}
#ifdef DDC
  , RaXi {RaXi_in}
  , tau {tau_in}
#endif
{
  calculateDerivedConstants();
}

void Constants::calculateDerivedConstants() {
  nX = nZ*aspectRatio;
  dz = 1.0/(nZ-1);
  dx = double(aspectRatio)/(nX-1);
  oodz2 = pow(1.0/dz, 2);
}

void Constants::print() const {
  std::cout <<"Parameters:" << std::endl;
  std::cout << "nZ: " << nZ << std::endl;
  std::cout << "nN: " << nN << std::endl;
  std::cout << "aspectRatio: " << aspectRatio << std::endl;
  std::cout << "Ra: " << Ra << std::endl;
#ifdef DDC
  std::cout << "RaXi: " << RaXi << std::endl;
  std::cout << "tau: " << tau << std::endl;
#endif
  std::cout << "Pr: " << Pr << std::endl;
  std::cout << "dt: " << initialDt << std::endl;
  std::cout << "totalTime: " << totalTime << std::endl;
  std::cout << "saveFolder: " << saveFolder << std::endl;
  std::cout << "icFile: " << icFile << std::endl;
}

bool Constants::isValid() const {
  if(nZ <=0 or nN <=0 or aspectRatio <= 0) {
    std::cout << " nZ (" << nZ
    << ") nN (" << nN
    << ") aspectRatio (" << aspectRatio
    << ") should be positive integers.\n" << std::endl;
    return -1;
  }
  if(initialDt <= 0.0f
  or Ra <= 0.0f
#ifdef DDC
  or RaXi <= 0.0f
  or tau <= 0.0f
#endif
  or Pr <= 0.0f
  or totalTime <= 0.0f
  or timeBetweenSaves <= 0.0f) {
    std::cout << " initial dt (" << initialDt
    << ") Ra (" << Ra
#ifdef DDC
    << ") RaXi (" << Ra
    << ") tau (" << Ra
#endif
    << ") Pr (" << Pr
    << ") total time (" << totalTime
    << ") time between saves (" << timeBetweenSaves
    << ") should be positive decimals.\n" << std::endl;
    return -1;
  }
  if(saveFolder == "" or icFile == "") {
    std::cout <<"Save folder and initial conditions file should be present.\n" << std::endl;
    return -1;
  }
  return 1;
}

void Constants::readJson(const std::string &filePath) {
  nlohmann::json j;

  std::ifstream in(filePath);
  if(in.is_open()) {
    in >> j;
    in.close();
  } else {
    std::cerr << "Could not read constants from JSON file " << filePath << std::endl;
  }

  nZ = j["nZ"];
  nN = j["nN"];
  initialDt = j["initialDt"];
  Ra = j["Ra"];
  Pr = j["Pr"];
  aspectRatio = j["aspectRatio"];
  timeBetweenSaves = j["timeBetweenSaves"];
  totalTime = j["totalTime"];
  saveFolder = j["saveFolder"];
  icFile = j["icFile"];
#ifdef DDC
  RaXi = j["RaXi"];
  tau = j["tau"];
#endif

  calculateDerivedConstants();
}

void Constants::writeJson(const std::string &filePath) const {
  nlohmann::json j;

  j["nZ"] = nZ;
  j["nN"] = nN;
  j["initialDt"] = initialDt;
  j["Ra"] = Ra;
  j["Pr"] = Pr;
  j["aspectRatio"] = aspectRatio;
  j["timeBetweenSaves"] = timeBetweenSaves;
  j["totalTime"] = totalTime;
  j["saveFolder"] = saveFolder;
  j["icFile"] = icFile;
#ifdef DDC
  j["RaXi"] = RaXi;
  j["tau"] = tau;
#endif

  std::ofstream out(filePath);
  if (out.is_open()) {
    out << j;
    out.close();
  } else {
    std::cerr << "Could not print constants to JSON file " << filePath << std::endl;
  }
}
