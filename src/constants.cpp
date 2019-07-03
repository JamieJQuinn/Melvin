#include <constants.hpp>

#include <json.hpp>
#include <fstream>
#include <iostream>

Constants::Constants() {}

Constants::~Constants() {}

Constants::Constants(const std::string &input) {
  readJson(input);
}

void Constants::calculateDerivedConstants() {
  nX = nZ*aspectRatio;
  dz = 1.0/(nZ-1);
  dx = real(aspectRatio)/(nX-1);
  oodz2 = pow(1.0/dz, 2);
}

void Constants::print() const {
  std::cout <<"Parameters:" << std::endl;
  std::cout << "nZ: " << nZ << std::endl;
  std::cout << "nN: " << nN << std::endl;
  std::cout << "aspectRatio: " << aspectRatio << std::endl;
  std::cout << "Ra: " << Ra << std::endl;
  std::cout << "Pr: " << Pr << std::endl;
  std::cout << "dt: " << initialDt << std::endl;
  std::cout << "totalTime: " << totalTime << std::endl;
  std::cout << "saveFolder: " << saveFolder << std::endl;
  std::cout << "is nonlinear? " << isNonlinear << std::endl;
  std::cout << "is double diffusion? " << isDoubleDiffusion << std::endl;
  std::cout << "is CUDA enabled? " << isCudaEnabled << std::endl;
  std::cout << "icFile: " << icFile << std::endl;
  std::cout << "CUDA threads per x: " << threadsPerBlock_x << std::endl;
  std::cout << "CUDA threads per y: " << threadsPerBlock_y << std::endl;

  if(isDoubleDiffusion) {
    std::cout << "RaXi: " << RaXi << std::endl;
    std::cout << "tau: " << tau << std::endl;
  }
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
  or Pr <= 0.0f
  or totalTime <= 0.0f
  or timeBetweenSaves <= 0.0f) {
    std::cout << " initial dt (" << initialDt
    << ") Ra (" << Ra
    << ") Pr (" << Pr
    << ") total time (" << totalTime
    << ") time between saves (" << timeBetweenSaves
    << ") should be positive decimals.\n" << std::endl;
    return -1;
  }

  if(isDoubleDiffusion) {
    if(RaXi <= 0.0f or tau <= 0.0f) {
      std::cout
      << ") RaXi (" << RaXi
      << ") tau (" << tau
      << ") should be positive floats.\n" << std::endl;
    }
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
  isNonlinear = j["isNonlinear"];
  isDoubleDiffusion = j["isDoubleDiffusion"];
  if (j.find("isCudaEnabled") != j.end()) {
    isCudaEnabled = j["isCudaEnabled"];
  } else {
    isCudaEnabled = false;
  }
  if(isCudaEnabled) {
    threadsPerBlock_x = j["threadsPerBlock_x"];
    threadsPerBlock_y = j["threadsPerBlock_y"];
  }
  icFile = j["icFile"];

  if(isDoubleDiffusion) {
    RaXi = j["RaXi"];
    tau = j["tau"];
  }

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
  j["isNonlinear"] = isNonlinear;
  j["isCudaEnabled"] = isCudaEnabled;
  j["icFile"] = icFile;
  if(isCudaEnabled) {
    j["threadsPerBlock_x"] = threadsPerBlock_x;
    j["threadsPerBlock_y"] = threadsPerBlock_y;
  }
  if(isDoubleDiffusion) {
    j["RaXi"] = RaXi;
    j["tau"] = tau;
  }

  std::ofstream out(filePath);
  if (out.is_open()) {
    out << j;
    out.close();
  } else {
    std::cerr << "Could not print constants to JSON file " << filePath << std::endl;
  }
}
