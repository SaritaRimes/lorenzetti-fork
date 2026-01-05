#include "CaloCluster/CaloCluster.h"
#include "G4Kernel/CaloPhiRange.h"
#include "ShowerShapes.h"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"
using namespace std;
using namespace SG;
ShowerShapes::ShowerShapes( std::string name ) : 
  IMsgService(name),
  AlgTool()
{;}
//!=====================================================================
StatusCode ShowerShapes::initialize()
{
  return StatusCode::SUCCESS;
}
//!=====================================================================
StatusCode ShowerShapes::finalize()
{
  return StatusCode::SUCCESS;
}
//!=====================================================================
StatusCode ShowerShapes::execute( SG::EventContext &/*ctx*/, Gaugi::EDM *edm ) const
{
  MSG_DEBUG("Calculate shower shapes for this cluster." );
  
  auto *clus = static_cast<xAOD::CaloCluster*>(edm);
  // Eratio for strip em layer (EM1, all cells must be EMB1 or EMEC1 sampling)
  auto em1Cells = clus->cells();
	em1Cells.erase(std::remove_if(em1Cells.begin(),em1Cells.end(),
                [](const xAOD::CaloCell* c){return (
                  (c->descriptor()->sampling() != CaloSampling::EMB1) &&
                  (c->descriptor()->sampling() != CaloSampling::EMEC1) );}),
                em1Cells.end());

  std::sort(em1Cells.begin(), em1Cells.end(),
            [](const xAOD::CaloCell* c1, const xAOD::CaloCell* c2){return c1->e() > c2->e();});
  
  
  float emaxs1  = (em1Cells.size()>=4)?(em1Cells[0]->e() + em1Cells[1]->e()):0;
  float e2tsts1 = (em1Cells.size()>=4)?(em1Cells[2]->e() + em1Cells[3]->e()):0;
  float eratio = (emaxs1 + e2tsts1)?((emaxs1 - e2tsts1)/(emaxs1 + e2tsts1)):0.;

  float e277  = sumEnergyEM( clus, 2, 7, 7 );
  float e233  = sumEnergyEM( clus, 2, 3, 3 );
  float e237  = sumEnergyEM( clus, 2, 3, 7 );
  float reta  = e237/e277;
  float rphi  = e233/e237;
  float e0    = sumEnergyEM( clus, 0 );
  float e1    = sumEnergyEM( clus, 1 );
  float e2    = sumEnergyEM( clus, 2 );
  float e3    = sumEnergyEM( clus, 3 );
  float ehad1 = sumEnergyHAD( clus, 0 );
  float ehad2 = sumEnergyHAD( clus, 1 );
  float ehad3 = sumEnergyHAD( clus, 2 );
  float weta2 = calculateWeta2( clus, 3, 5 );
  // distance in eta and phi between the energy weighted cluster centroid and the hottest cell
  float detaCentHotCell = calculateDetaEnergeticCentroidHotCell( clus );
  float dphiCentHotCell = calculateDphiEnergeticCentroidHotCell( clus );
  // distance in eta and phi between the energy weighted cluster centroid and the cell where it is located
  float detaCentCell = calculateDetaEnergeticCentroidCell( clus );
  float dphiCentCell = calculateDphiEnergeticCentroidCell( clus );

  float tofClus = calculateClusterTof( clus );

  float etot = e0+e1+e2+e3+ehad1+ehad2+ehad3;
  float emtot = e0+e1+e2+e3;
  // fraction of energy deposited in 1st sampling
  float f0 = e0 / emtot;
  float f1 = e1 / emtot;
  float f2 = e2 / emtot;
  float f3 = e3 / emtot;
  float rhad = (ehad1+ehad2+ehad3) / emtot;
  float rhad1 = (ehad1) / emtot;

  clus->setE0( e0 );
  clus->setE1( e1 );
  clus->setE2( e2 );
  clus->setE3( e3 );
  clus->setEhad1( ehad1 );
  clus->setEhad2( ehad2 );
  clus->setEhad3( ehad3 );
  clus->setEtot( etot );
  clus->setE277( e277 );
  clus->setE237( e237 );
  clus->setE233( e233 );
  clus->setReta( reta );
  clus->setRphi( rphi );
  clus->setWeta2( weta2 );
  clus->setDetaCentHotCell( detaCentHotCell );
  clus->setDphiCentHotCell( dphiCentHotCell );
  clus->setDetaCentCell( detaCentCell );
  clus->setDphiCentCell( dphiCentCell );
  clus->setTofClus( tofClus );
  // Eratio for strip em layer (EM1)
  clus->setEratio( eratio );
  // Extra eratio information for each sampling layer
  clus->setEmaxs1( emaxs1 );
  clus->setE2tsts1( e2tsts1 );
  clus->setF0( f0 );
  clus->setF1( f1 );
  clus->setF2( f2 );
  clus->setF3( f3 );
  clus->setRhad( rhad );
  clus->setRhad1( rhad1 );
  // Only EM energy since this is a eletromagnetic cluster
  clus->setEt( clus->eta() != 0.0 ? (e0+e1+e2+e3)/cosh(std::abs(clus->eta())) : 0.0 ); 
  clus->setE(e0+e1+e2+e3);


  if(clus->isForward()){
    std::vector<TVector3> axis = calculateShowerAxis(clus);
    float secondR = calculateSecondR(clus, axis);
    float lambdaCenter = calculateLambdaCenter(clus, axis);
    float secondLambda = calculateSecondLambda(clus, axis);
    float fracMax = calculateFracMax(clus, axis);
    float lateralMom = calculateLateralMom(clus, axis);
    float longitudinalMom = calculateLongitudinalMom(clus, axis);
    clus->setSecondR(secondR);
    clus->setLambdaCenter(lambdaCenter);
    clus->setSecondLambda(secondLambda);
    clus->setFracMax(fracMax);
    clus->setLateralMom(lateralMom);
    clus->setLongitudinalMom(longitudinalMom);
  }
  
  return StatusCode::SUCCESS;
}

//!=====================================================================
float ShowerShapes::sumEnergyEM( xAOD::CaloCluster *clus, int sampling, unsigned eta_ncell, unsigned phi_ncell ) const
{
  float energy = 0.0;
  for ( const auto& cell : clus->cells() )
  {
    const xAOD::CaloDetDescriptor* det = cell->descriptor();
    // Exclude all HEC, TILE cells from the loop...
    if (det->detector()!= Detector::TTEM) continue;
    // if PS, the cell must be PSE or PSB
    if(sampling==0 && (
      det->sampling()!=CaloSampling::PSB &&
      det->sampling()!=CaloSampling::PSE )  ) continue;
    // if EM1, this cell must be EMB1 or EMEC1
    if(sampling==1 && (
      det->sampling()!=CaloSampling::EMB1 &&
      det->sampling()!=CaloSampling::EMEC1 )  ) continue;
    // if EM2, this cell must be EMB2 or EMEC2
    if(sampling==2 && (
      det->sampling()!=CaloSampling::EMB2 &&
      det->sampling()!=CaloSampling::EMEC2 )  ) continue;
    // if EM3, this cell must be EMB3 or EMEC3
    if(sampling==3 && (
      det->sampling()!=CaloSampling::EMB3 &&
      det->sampling()!=CaloSampling::EMEC3 )  ) continue;
    float deltaEta = std::abs( clus->eta() - cell->eta() );
    float deltaPhi = std::abs( CaloPhiRange::fix( clus->phi() - cell->phi() ) );
    
    if( deltaEta < eta_ncell*(cell->deltaEta()/2) && deltaPhi < phi_ncell*(cell->deltaPhi()/2) ){
      energy+= cell->e();
    }
  }
  return energy;
}
//!=====================================================================
float ShowerShapes::sumEnergyHAD( xAOD::CaloCluster *clus, int sampling, unsigned eta_ncell, unsigned phi_ncell ) const
{
  float energy = 0.0;
  for ( const auto& cell : clus->cells() )
  {
    const xAOD::CaloDetDescriptor* det = cell->descriptor();
    // Exclude all EM cells from the loop...
    if (det->detector()== Detector::TTEM) continue;
    // HAD 1
    if(sampling==0 && 
    ( (det->sampling()!=CaloSampling::HEC1) &&
      (det->sampling()!=CaloSampling::TileCal1) &&
      (det->sampling()!=CaloSampling::TileExt1) ) ) continue;
    // HAD 2
    if(sampling==1 && 
    ( (det->sampling()!=CaloSampling::HEC2) &&
      (det->sampling()!=CaloSampling::TileCal2) &&
      (det->sampling()!=CaloSampling::TileExt2) ) ) continue;
    // HAD 3
    if(sampling==2 && 
    ( (det->sampling()!=CaloSampling::HEC3) &&
      (det->sampling()!=CaloSampling::TileCal3) &&
      (det->sampling()!=CaloSampling::TileExt3) ) ) continue;
    float deltaEta = std::abs( clus->eta() - cell->eta() );
    float deltaPhi = std::abs( CaloPhiRange::fix( clus->phi() - cell->phi() ) );
    
    if( deltaEta < eta_ncell*(cell->deltaEta()/2) && deltaPhi < phi_ncell*(cell->deltaPhi()/2) ){
      energy+= cell->e();
    }
  }
  return energy;
}
//!=====================================================================
float ShowerShapes::calculateWeta2( xAOD::CaloCluster *clus , unsigned eta_ncell, unsigned phi_ncell) const
{
  float En2=0.0;
  float En=0.0;
  float E=0.0;
  for ( auto& cell : clus->cells() ){
    const xAOD::CaloDetDescriptor* det = cell->descriptor();
    // Only EM2 cells should be used here...
    if(det->sampling()!=CaloSampling::EMB2 && det->sampling()!=CaloSampling::EMEC2)  continue;
    
    float deltaEta = std::abs( clus->eta() - cell->eta() );
    float deltaPhi = std::abs( CaloPhiRange::diff( clus->phi() , cell->phi() ) );
 
    if( deltaEta < eta_ncell*(cell->deltaEta()/2) && deltaPhi < phi_ncell*(cell->deltaPhi()/2) ){
      En2 += cell->e() * std::pow(cell->eta(),2);
      En += cell->e() * cell->eta();
      E += cell->e();
    }
  }// Loop over all cells inside of the cluster
  return std::sqrt( (En2/E) - std::pow( (En/E),2 ) );
}

std::vector<TVector3> ShowerShapes::calculateShowerAxis(xAOD::CaloCluster * clus) const{
  float xc = 0;
  float yc = 0;
  float zc = 0;
  float w = 0;
  float emax1 = 0;
  float emax2 = 0;
  double angle = 0;
  float totalE = 0;
  float m_maxAxisAngle = 20*3.14/180;
  for ( auto& cell : clus->cells() ){
    float eta = cell->eta();
    float phi = cell->phi();
    float energy = cell->e();
	  // gets the 2 highest cell energies in the cluster
	  if (energy > emax2){
		  emax2 = energy;
    }
		if (energy > emax1){
			emax2 = emax1;
			emax1 = energy;
		}

    float thetaCell = 2*atan(pow(2.718281828, -eta));
    float zCell = cell->descriptor()->z();
    float rCell = zCell/cos(thetaCell);
    float xCell = rCell*sin(thetaCell)*cos(phi);
    float yCell = rCell*sin(thetaCell)*sin(phi);
    xc += energy*xCell;
    yc += energy*yCell;
    zc += energy*zCell;
    totalE+=energy;
  }
  TVector3 emec1_cell(0,0,0);
  TVector3 showerAxis(xc/totalE,yc/totalE,zc/totalE);
  TVector3 showerCenter(xc/totalE,yc/totalE,zc/totalE);
  showerAxis.SetMag(1.0);
  TMatrixD C(3,3);
  for ( auto& cell : clus->cells() ){
    float energy = cell->e();
    float eta = cell->eta();
    float phi = cell->phi();
    float thetaCell = 2*atan(pow(2.718281828, -eta));
    const double e2 = energy * energy;
    float zCell = cell->descriptor()->z();
    float rCell = zCell/cos(thetaCell);
    float xCell = rCell*sin(thetaCell)*cos(phi);
    float yCell = rCell*sin(thetaCell)*sin(phi);
    C(0,0) += e2*(xCell-xc)*(xCell-xc);
    C(1,0) += e2*(xCell-xc)*(yCell-yc);
    C(2,0) += e2*(xCell-xc)*(zCell-zc);
    C(1,1) += e2*(yCell-yc)*(yCell-yc);
    C(2,1) += e2*(yCell-yc)*(zCell-zc);
    C(2,2) += e2*(zCell-zc)*(zCell-zc);
    w += e2;
	
	  if (cell->descriptor()->sampling() == CaloSampling::EMEC1){
		  emec1_cell.SetZ(zCell);
	  }
  }

  C(0,0) = C(0,0) /w;
  C(1,0) = C(1,0) /w;
  C(2,0) = C(2,0) /w;
  C(1,1) = C(1,1) /w;
  C(2,1) = C(2,1) /w;
  C(2,2) = C(2,2) /w;
  auto S = TMatrixDEigen(C).GetEigenValues() ;
  auto U = TMatrixDEigen(C).GetEigenVectors() ;
  float epsilon = 1.E-6;
  if ( std::abs(S(0,0)) >= epsilon && std::abs(S(1,1)) >= epsilon && std::abs(S(2,2)) >= epsilon ) { 
    TVector3 prAxis(showerAxis);
    int iEigen = -1;
    for (int i=0;i<3;i++) {
      TVector3 tmpAxis(U(i,0), U(i,1), U(i,2));
      double tmpAngle = tmpAxis.Angle(showerAxis);
      if ( tmpAngle > 90*deg ) { 
        tmpAngle = 180*deg - tmpAngle;
        tmpAxis = -tmpAxis;
      }
       if ( iEigen == -1 || tmpAngle < angle ) {
        iEigen = i;
        angle = tmpAngle;
        prAxis = tmpAxis;
       }
      if ( angle < m_maxAxisAngle ) {
        showerAxis = prAxis;
      }
    }
  }
	TVector3 energies(totalE, emax1, emax2);
  std::vector<TVector3> axis = {showerAxis, showerCenter, energies, emec1_cell};
  return axis;
}

float ShowerShapes::calculateLambdaCenter ( xAOD::CaloCluster *clus, std::vector<TVector3> axis ) const{
  float LambdaCenter = 0;
  float factor = 0;
	auto showerAxis = axis.at(0);
	auto z_showerAxis = axis.at(0).Z();
	auto showerCenter = axis.at(1);
	auto z_showerCenter = axis.at(1).Z();
	auto z_emec1 = axis.at(3).Z();

	if (z_showerCenter > 0) {
		factor = (z_showerCenter - z_emec1)/z_showerAxis;
	}
	else {
		factor = (z_showerCenter + z_emec1)/z_showerAxis;
	}
	LambdaCenter = (factor*showerAxis).Mag();
	
	return LambdaCenter;
}

float ShowerShapes::calculateSecondLambda ( xAOD::CaloCluster *clus, std::vector<TVector3> axis ) const{
  float secondLambda = 0;
	auto showerAxis = axis.at(0);
	auto showerCenter = axis.at(1);
	auto totalE = axis.at(2).X();
    
	for ( auto& cell : clus->cells() ){
    	float eta = cell->eta();
      	float phi = cell->phi();
      	float energy = cell->e();
      	float thetaCell = 2*atan(pow(2.718281828, -eta));
      	float zCell = cell->descriptor()->z();
      	float rCell = zCell/cos(thetaCell);
      	float xCell = rCell*sin(thetaCell)*cos(phi);
      	float yCell = rCell*sin(thetaCell)*sin(phi);
      	TVector3 currentCell(xCell,yCell,zCell);
      	float r = ((currentCell-showerCenter).Dot(showerAxis));
      	// if (cell->e() < 0) continue;
      	secondLambda+=energy*r*r;
  	}
	
	return secondLambda/totalE;
}

float ShowerShapes::calculateSecondR ( xAOD::CaloCluster *clus, std::vector<TVector3> axis ) const{
  float secondR = 0;
	auto showerAxis = axis.at(0);
	auto showerCenter = axis.at(1);
	auto totalE = axis.at(2).X();
    
	for ( auto& cell : clus->cells() ){
    	  float eta = cell->eta();
      	float phi = cell->phi();
      	float energy = cell->e();
      	float thetaCell = 2*atan(pow(2.718281828, -eta));
      	float zCell = cell->descriptor()->z();
      	float rCell = zCell/cos(thetaCell);
      	float xCell = rCell*sin(thetaCell)*cos(phi);
      	float yCell = rCell*sin(thetaCell)*sin(phi);
      	TVector3 currentCell(xCell,yCell,zCell);
      	float r = ((currentCell-showerCenter).Cross(showerAxis)).Mag();
      	// if (cell->e() < 0) continue;
      	secondR+=energy*r*r;
  	}	
	return secondR/totalE;
}

float ShowerShapes::calculateFracMax ( xAOD::CaloCluster *clus, std::vector<TVector3> axis ) const{
  float fracMax = 0;
	auto totalE = axis.at(2).X();
	auto emax1 = axis.at(2).Y();
	fracMax = emax1;
	return fracMax/totalE;
}

float ShowerShapes::calculateLateralMom ( xAOD::CaloCluster *clus, std::vector<TVector3> axis ) const{
	float lateralMom = 0;
  float r = 0;
  float r_2 = 0;
  float r_max = 0;
	auto showerAxis = axis.at(0);
	auto showerCenter = axis.at(1);
	auto totalE = axis.at(2).X();
	auto emax1 = axis.at(2).Y();
	auto emax2 = axis.at(2).Z();
  
  for ( auto& cell : clus->cells() ){
      	float eta = cell->eta();
      	float phi = cell->phi();
      	float energy = cell->e();
      	float thetaCell = 2*atan(pow(2.718281828, -eta));
      	float zCell = cell->descriptor()->z();
      	float rCell = zCell/cos(thetaCell);
      	float xCell = rCell*sin(thetaCell)*cos(phi);
      	float yCell = rCell*sin(thetaCell)*sin(phi);
      	TVector3 currentCell(xCell,yCell,zCell);
	  	if (energy != emax1 && energy != emax2){
	  		r = ((currentCell-showerCenter).Cross(showerAxis)).Mag();
			  r_2 += energy*pow(r,2)/totalE;
	  	}
	  	else {
			r_max += energy*pow(40,2)/totalE;
      }
	}
  lateralMom = r_2/(r_2+r_max);
  return lateralMom;
}

float ShowerShapes::calculateLongitudinalMom ( xAOD::CaloCluster *clus, std::vector<TVector3> axis ) const{
	float longitudinalMom = 0;
  float l = 0;
  float l_2 = 0;
  float l_max = 0;
	auto showerAxis = axis.at(0);
	auto showerCenter = axis.at(1);
	auto totalE = axis.at(2).X();
	auto emax1 = axis.at(2).Y();
	auto emax2 = axis.at(2).Z();
  
    for ( auto& cell : clus->cells() ){
      float eta = cell->eta();
      float phi = cell->phi();
      float energy = cell->e();
      float thetaCell = 2*atan(pow(2.718281828, -eta));
      float zCell = cell->descriptor()->z();
      float rCell = zCell/cos(thetaCell);
      float xCell = rCell*sin(thetaCell)*cos(phi);
      float yCell = rCell*sin(thetaCell)*sin(phi);
      TVector3 currentCell(xCell,yCell,zCell);
	  	if (energy != emax1 && energy != emax2){
	  		l = ((currentCell-showerCenter).Dot(showerAxis));
			  l_2 += energy*pow(l,2)/totalE;
	  	}
	  	else {
			  l_max += energy*pow(100,2)/totalE;
	  	} 
  	}
    longitudinalMom = l_2/(l_2+l_max);
  	return longitudinalMom;
}
//!=====================================================================
float ShowerShapes::calculateDetaEnergeticCentroidHotCell( xAOD::CaloCluster *clus ) const{
  float sumEtaEnergy = 0, sumEnergy = 0;
  float etaMostEnergCell = 0, eMostEnergCell = -1;

  // Estimating the eta coordinate of the energetic centroid
  for(auto& cell : clus->cells()) {
    float energy = cell->e();

    if(energy < 50) continue;

    float eta = cell->eta();

    sumEtaEnergy += eta * energy;
    sumEnergy += energy;

    if(energy > eMostEnergCell) {
      etaMostEnergCell = eta;
      eMostEnergCell = energy;
    }    
  }

  if(sumEnergy == 0) {
    std::cerr << "Warning: sum of energy is zero." << std::endl;
  }

  float etaCent = sumEtaEnergy/sumEnergy;
  float detaCentHotCell = etaCent - etaMostEnergCell;
  
  return detaCentHotCell;
}
//!=====================================================================
float ShowerShapes::calculateDetaEnergeticCentroidCell( xAOD::CaloCluster *clus ) const{
  float sumEtaEnergy = 0, sumEnergy = 0;
  
  // Estimating the eta coordinate of the energetic centroid
  for(auto& cell : clus->cells()) {
    float energy = cell->e();

    if(energy < 50) continue;

    float eta = cell->eta();

    sumEtaEnergy += eta * energy;
    sumEnergy += energy;
  }

  if(sumEnergy == 0) {
    std::cerr << "Warning: sum of energy is zero." << std::endl;
  }

  float etaCent = sumEtaEnergy/sumEnergy;

  // Finding the cell where the centroid is located
  float etaCellCentroid;
  for(auto& cell : clus->cells()) {
    float eta = cell->eta();
    float deta = cell->deltaEta();

    if(std::fabs(etaCent - eta) <= deta/2) {
      etaCellCentroid = eta;
      break;
    }
  }

  float detaCentCell = etaCent - etaCellCentroid;
  
  return detaCentCell;
}
//!=====================================================================
float ShowerShapes::calculateDphiEnergeticCentroidHotCell( xAOD::CaloCluster *clus ) const{
  const float pi = 3.141593, threshold = 5*0.025 /* approx pi/25 */; 
  
  // Checking if there are cells in the phi transition region (-pi, +pi)
  float phiMin = 3.5; // higher than pi
  float phiMax = -3.5; // smaller than -pi
  for(auto& cell : clus->cells()) {
    float phi = CaloPhiRange::fix(cell->phi());
    
    if(phi < phiMin) { phiMin = phi; }
    if(phi > phiMax) { phiMax = phi; }
  }
  
  bool inTransitionRegion = std::fabs(phiMax - phiMin) > (2*pi - threshold);

  // Estimating the phi coordinate of the energetic centroid
  float phiMostEnergCell = 0, eMostEnergCell = -1;
  float phiCent; 
  if(inTransitionRegion) { // if it is the transition region
    float sumPhiEnergyPos = 0, sumPhiEnergyNeg = 0;
    float sumEnergyPos = 0, sumEnergyNeg = 0;

    for(auto& cell : clus->cells()) {
      float energy = cell->e();

      if(energy < 50) continue;

      float phi = CaloPhiRange::fix(cell->phi());

      if(phi >= 0) {
        sumPhiEnergyPos += phi * energy;
        sumEnergyPos += energy;
      } else {
        sumPhiEnergyNeg += phi * energy; 
        sumEnergyNeg += energy;
      }

      if(energy > eMostEnergCell) {
        phiMostEnergCell = phi;
        eMostEnergCell = energy;
      }
    }

    if(sumEnergyPos == 0) {
      std::cerr << "Warning: sum of energy in the positive side of phi is zero." << std::endl;
    }

    if(sumEnergyNeg == 0) {
      std::cerr << "Warning: sum of energy in the negativa side of phi is zero." << std::endl;
    }

    // Estimating the weighted centroid for both positive and negative part
    float phiCentPos = CaloPhiRange::fix(sumPhiEnergyPos/sumEnergyPos);
    float phiCentNeg = CaloPhiRange::fix(sumPhiEnergyNeg/sumEnergyNeg);

    // Calculating the phi of the centroid using the barycenter and considering phi_neg was 0
    float epsilonPos = CaloPhiRange::fix(pi - phiCentPos);
    float epsilonNeg = CaloPhiRange::fix(phiCentNeg + pi);
    phiCent = CaloPhiRange::fix(((epsilonPos + epsilonNeg) * (sumEnergyPos/(sumEnergyPos + sumEnergyNeg))) + phiCentNeg);

  } else { // If it is not the transition region
    float sumPhiEnergy = 0, sumEnergy = 0;

    for(auto& cell : clus->cells()) {
      float energy = cell->e();

      if(energy < 50) continue;

      float phi = CaloPhiRange::fix(cell->phi());

      sumPhiEnergy += phi * energy;
      sumEnergy += energy;

      if(energy > eMostEnergCell) {
        phiMostEnergCell = phi;
        eMostEnergCell = energy;
      }
    }

    if(sumEnergy == 0) {
      std::cerr << "Warning: sum of energy is zero." << std::endl;
    }

    phiCent = CaloPhiRange::fix(sumPhiEnergy/sumEnergy);
  }
  
  float dphiCentHotCell = CaloPhiRange::diff(phiCent, phiMostEnergCell);

  return dphiCentHotCell;
}
// !=====================================================================
float ShowerShapes::calculateDphiEnergeticCentroidCell( xAOD::CaloCluster *clus ) const{
  const float pi = 3.141593, threshold = 5*0.025 /* approx pi/25 */, dphiTile2 = 0.025; 
  
  // Checking if there are cells in the phi transition region (-pi, +pi)
  float phiMin = 3.5; // higher than pi
  float phiMax = -3.5; // smaller than -pi
  for(auto& cell : clus->cells()) {
    float phi = CaloPhiRange::fix(cell->phi());
    
    if(phi < phiMin) { phiMin = phi; }
    if(phi > phiMax) { phiMax = phi; }
  }
  
  bool inTransitionRegion = std::fabs(phiMax - phiMin) > (2*pi - threshold);

  // Estimating the phi coordinate of the energetic centroid
  float phiCent; 
  if(inTransitionRegion) { // If it is the transition region
    float sumPhiEnergyPos = 0, sumPhiEnergyNeg = 0;
    float sumEnergyPos = 0, sumEnergyNeg = 0;

    for(auto& cell : clus->cells()) {
      float energy = cell->e();

      if(energy < 50) continue;

      float phi = CaloPhiRange::fix(cell->phi());

      if(phi >= 0) {
        sumPhiEnergyPos += phi * energy;
        sumEnergyPos += energy;
      } else {
        sumPhiEnergyNeg += phi * energy; 
        sumEnergyNeg += energy;
      }
    }

    if(sumEnergyPos == 0) {
      std::cerr << "Warning: sum of energy in the positive side of phi is zero." << std::endl;
    }

    if(sumEnergyNeg == 0) {
      std::cerr << "Warning: sum of energy in the negativa side of phi is zero." << std::endl;
    }

    // Estimating the weighted centroid for both positive and negative part
    float phiCentPos = CaloPhiRange::fix(sumPhiEnergyPos/sumEnergyPos);
    float phiCentNeg = CaloPhiRange::fix(sumPhiEnergyNeg/sumEnergyNeg);

    // Calculating the phi of the centroid using the barycenter and considering phi_neg was 0
    float epsilonPos = CaloPhiRange::fix(pi - phiCentPos);
    float epsilonNeg = CaloPhiRange::fix(phiCentNeg + pi);
    phiCent = CaloPhiRange::fix(((epsilonPos + epsilonNeg) * (sumEnergyPos/(sumEnergyPos + sumEnergyNeg))) + phiCentNeg);

  } else { // If it is not the transition region
    float sumPhiEnergy = 0, sumEnergy = 0;

    for(auto& cell : clus->cells()) {
      float energy = cell->e();

      if(energy < 50) continue;

      float phi = CaloPhiRange::fix(cell->phi());

      sumPhiEnergy += phi * energy;
      sumEnergy += energy;
    }

    if(sumEnergy == 0) {
      std::cerr << "Warning: sum of energy is zero." << std::endl;
    }

    phiCent = CaloPhiRange::fix(sumPhiEnergy/sumEnergy);
  }

  // Finding the cell where the centroid is located
  float phiCellCentroid;
  int numberCellCentroid;
  if(phiCent >= 0) {
    numberCellCentroid = (int)((pi - phiCent)/dphiTile2);
    phiCellCentroid = pi - (numberCellCentroid * dphiTile2) - (dphiTile2/2);
  } else {
    numberCellCentroid = (int)((phiCent - (-pi))/dphiTile2);
    phiCellCentroid = -pi + (numberCellCentroid * dphiTile2) + (dphiTile2/2); 
  }
  
  float dphiCentCell = CaloPhiRange::diff(phiCent, phiCellCentroid);

  return dphiCentCell;
}
// !=====================================================================
float ShowerShapes::calculateClusterTof( xAOD::CaloCluster *clus ) const {
  std::ofstream file("test_tof.csv");

  for(auto& cell : clus->cells()) {
    std::cout << "\n\n----- Inicia celula -----" << std::endl;
    
    // retrieving descriptor information from the cell
    const xAOD::CaloDetDescriptor* descri = cell->descriptor();

    int bcid_start = descri->bcid_start();
    int bcid_end = descri->bcid_end();

    for(int i = bcid_start; i <= bcid_end ; i++){
      std::cout << i << ", ";
    }
    std::cout << std::endl;

    file << descri->hash() << ", ";
    for(int i = bcid_start; i <= bcid_end ; i++){
      std::cout << "tof(" << i << "): " << descri->tof(i) << ", ";
      std::cout << "hash: " << descri->hash() << ", ";
      std::cout << "eta: " << descri->eta() << ", ";
      std::cout << "phi: " << descri->phi() << std::endl;

      file << descri->tof(i);
      if(i != bcid_end) file << ", ";
    }
    file << "\n";

    // tentar salvar dados de tof no csv pra plotar

    // auto tof = cell->descriptor()->tof(-8);


    // std::cout << "\ntof: " << tof << std::endl;
    // Encontrar link entre a celula do objeto CaloCell e a do objeto CaloDetDescriptor
    // Preciso identificar quais células do Descriptor estão no cluster de entrada
    // O CaloDetDescriptor tem um hash e o CaloCell tem um descriptor_link, a partir deles da pra mapear 
      // as celulas, pois os numeros sao os mesmos nos dois, ligando a celula em questao
      
    std::cout << "----- Finaliza celula -----" << std::endl;
  }
  file.close();
}