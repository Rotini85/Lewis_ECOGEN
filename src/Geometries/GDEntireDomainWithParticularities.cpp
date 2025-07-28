//  
//       ,---.     ,--,    .---.     ,--,    ,---.    .-. .-. 
//       | .-'   .' .')   / .-. )  .' .'     | .-'    |  \| | 
//       | `-.   |  |(_)  | | |(_) |  |  __  | `-.    |   | | 
//       | .-'   \  \     | | | |  \  \ ( _) | .-'    | |\  | 
//       |  `--.  \  `-.  \ `-' /   \  `-) ) |  `--.  | | |)| 
//       /( __.'   \____\  )---'    )\____/  /( __.'  /(  (_) 
//      (__)              (_)      (__)     (__)     (__)     
//      Official webSite: https://code-mphi.github.io/ECOGEN/
//
//  This file is part of ECOGEN.
//
//  ECOGEN is the legal property of its developers, whose names 
//  are listed in the copyright file included with this source 
//  distribution.
//
//  ECOGEN is free software: you can redistribute it and/or modify
//  it under the terms of the GNU General Public License as published 
//  by the Free Software Foundation, either version 3 of the License, 
//  or (at your option) any later version.
//  
//  ECOGEN is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
//  GNU General Public License for more details.
//  
//  You should have received a copy of the GNU General Public License
//  along with ECOGEN (file LICENSE).  
//  If not, see <http://www.gnu.org/licenses/>.

#include "GDEntireDomainWithParticularities.h"

//***********************************************************


GDEntireDomainWithParticularities::GDEntireDomainWithParticularities(std::string name, std::vector<Phase*> vecPhases, Mixture* mixture, std::vector<Transport> vecTransports, const int& physicalEntity, const std::string& icType)
  : GeometricalDomain(name, vecPhases, mixture, vecTransports, physicalEntity), m_icType(icType)
{
  if (m_icType == "ic_Rayleigh-TaylorInstability") {
    initializeRayleighTaylorInterface();
  }
}
void GDEntireDomainWithParticularities::initializeRayleighTaylorInterface() {
  const double pi = std::acos(-1);
  double lambda = 0.2;
  double h = 0.7;
  double k = 2 * pi / lambda;
  int nx = 1000;
  double amp = 0.05 / k;
  m_interfaceX.clear();
  m_interfaceY.clear();
  m_interfaceX.reserve(nx);
  m_interfaceY.reserve(nx);
  for (int i = 0; i < nx; i++) {
    double x = i * lambda / (nx - 1.);
    m_interfaceX.push_back(x);
    m_interfaceY.push_back(amp * std::sin(2. * pi * x / lambda + pi / 2.) + h);
  }
}

//***********************************************************

GDEntireDomainWithParticularities::~GDEntireDomainWithParticularities() {}

//***********************************************************

bool GDEntireDomainWithParticularities::belong(Coord& /*posElement*/, const int& /*lvl*/) const
{
  //1. Laplace pressure initialization
  //----------------------------------
  return true; //always belong to entire domain

  //2. Respecting special coordinates
  //---------------------------------
  //bool result(false);
  //if (posElement.getY() - posElement.getX() >= -2.-1.e-8) { result = true; }
  //return result;

  //3. Respecting special coordinates (with AMR test)
  //-------------------------------------------------
  // bool result(false);
  // if (lvl > 0) {
  //   if (posElement.getX() < 0.02 / std::pow(2., (double)(lvl))) { result = true; }
  // }
  // else {
  //   if (posElement.getX() < 0.02) { result = true; }
  // }
  // return result;

  //4. Random velocity perturbations
  //--------------------------------
  // return true; //always belong to entire domain

  //5. Rayleigh-Taylor instability
  //------------------------------
  // return true; //always belong to entire domain

  //6. Blast-wave equation
  //----------------------
  // return true; //always belong to entire domain
}

//******************************************************************

void GDEntireDomainWithParticularities::fillIn(Cell* cell) const
{
  //As basic fillIn: Test if the cell belongs to the geometrical domain
  bool belongs(true);
  if (cell->getElement() != 0) {
    Coord coord(cell->getPosition());
    if (!this->belong(coord, cell->getLvl())) { belongs = false; }
    //Test if the cell belongs to physical mesh entity (for unstructured meshes)
    if (cell->getElement()->getAppartenancePhysique() > 0 && m_physicalEntity > 0) {
      if (cell->getElement()->getAppartenancePhysique() != m_physicalEntity) { belongs = false; }
    }
  }

  if (belongs) {
    for (int k = 0; k < numberPhases; k++) { cell->copyPhase(k, m_vecPhases[k]); }
    cell->copyMixture(m_mixture);
    for (int k = 0; k < numberTransports; k++) { cell->setTransport(m_vecTransports[k].getValue(), k); }
    if(m_physicalEntity == -1){ cell->setWall(true); }
    else{ cell->setWall(false); }

    // Select IC logic based on m_icType
    // 1. Laplace pressure initialization
    if (m_icType == "ic_LaplacePressureInitialization") {
      if (cell->getElement() != 0) {
        double pressure(0.);
        Coord posElement(cell->getPosition());
        double radius;
        // 2D version
        radius = std::sqrt(std::pow(posElement.getX() - 2.e-4, 2.) + std::pow(posElement.getY(), 2.));
        pressure = 50.6625e5 + 1.e-4 / radius * (3.55e3 - 50.6625e5);
        for (int k = 0; k < numberPhases; k++) cell->getPhase(k)->setPressure(pressure);
        cell->getMixture()->setPressure(pressure);
      }
    }
    // 2. Respecting special coordinates
    else if (m_icType == "ic_Rayleigh-TaylorInstability") {
      if (cell->getElement() != 0) {
        Coord posElement(cell->getPosition());
        int index = 0;
        double minVal = 1e9;
        for (unsigned int i = 0; i < m_interfaceX.size(); i++) {
          if (std::abs(posElement.getX() - m_interfaceX[i]) < minVal) {
            minVal = std::abs(posElement.getX() - m_interfaceX[i]);
            index = i;
          }
        }
        if (index < 0 || index >= (int)m_interfaceY.size()) index = 0;
        double yInterface = m_interfaceY[index];
        double posY = posElement.getY();
        double rhoH = cell->getPhase(0)->getDensity();
        double rhoL = cell->getPhase(1)->getDensity();
        double pref = 1.e5, pinterface = pref, pressure = 0.;
        double g = 9.81, ly = 1.2;
        std::vector<double> alpha(2, 0.);
        if (posY > yInterface) {
          alpha[0] = 1.;
          alpha[1] = 0.;
          pressure = pref + rhoH * g * (ly - posY);
        } else {
          alpha[0] = 0.;
          alpha[1] = 1.;
          pinterface = pref + rhoH * g * (ly - yInterface);
          pressure = pinterface + rhoL * g * (yInterface - posY);
        }
        for (int k = 0; k < numberPhases; k++) {
          cell->getPhase(k)->setAlpha(alpha[k]);
          cell->getPhase(k)->setPressure(pressure);
        }
        cell->getMixture()->setPressure(pressure);
      }
    }

        int index = 0;
        
        if (cell->getElement() != 0) {
          Coord posElement(cell->getPosition());
          double minVal = 1.;
          for (unsigned int i = 0; i < interfaceX.size(); i++) { // Find nearest index corresponding to x-position of interface fn
            if (std::abs(posElement.getX() - interfaceX[i]) < minVal) { 
              minVal = std::abs(posElement.getX() - interfaceX[i]);
              index = i;
            }
          }

          // Check location to interface and initialize domain accordingly
          if (posElement.getY() > interfaceY[index]) {
          alpha[0] = 1.;
          alpha[1] = 0.;
          pressure = pref + rhoHeavy * g * (ly - posElement.getY());
          }
          else {
            alpha[0] = 0.;
            alpha[1] = 1.;
            pinterface = pref + rhoHeavy * g * (ly - interfaceY[index]);
            pressure = pinterface + rhoLight * g * (interfaceY[index] - posElement.getY());
          }

          for (int k = 0; k < numberPhases; k++) { 
            cell->getPhase(k)->setAlpha(alpha[k]);
            cell->getPhase(k)->setPressure(pressure);
          }
          cell->getMixture()->setPressure(pressure);
        }
      }
    }
    // 6. Blast-wave equation
    else if (m_icType == "ic_Blast-WaveEquation") {
      if (cell->getElement() != 0) {
        double pressure(0.), velocity(0.), pk(0.);
        double beta(1.48e6), omega(1.21e6);
        double posX(cell->getPosition().getX());
        double p0(1.01325e5), pS(35e6), soundSpeed(1625.), density(1000.);
        double shockFront(7.5e-3);
        if (posX < shockFront) {
          double r = posX - shockFront;
          pressure = p0 + 2 * pS * exp(beta * r / soundSpeed) * std::cos(- omega * r / soundSpeed + M_PI / 3.);
          velocity = (pressure - p0) / (density * soundSpeed);
          for (int k = 0; k < numberPhases; k++) {
            pk = pressure;
            cell->getPhase(k)->getEos()->verifyAndModifyPressure(pk);
            cell->getPhase(k)->setPressure(pk);
          }
          cell->getMixture()->setPressure(pressure);
          cell->getMixture()->setU(velocity);
        }
      }
    }
    // Default: do nothing or fallback logic
    else {
      // Default: do nothing or fallback logic
    }
  }
}
//***********************************************************
