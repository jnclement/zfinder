#include "zfinder.h"
#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>
#include <jetbase/JetContainerv1.h>
#include <jetbase/Jet.h>
#include <globalvertex/GlobalVertex.h>     // for GlobalVertex, GlobalVe...
#include <globalvertex/GlobalVertexMap.h> // for GlobalVertexMap
#include <globalvertex/GlobalVertexMapv1.h>
#include <globalvertex/GlobalVertexv2.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/RawTowerGeomContainer_Cylinderv1.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfoContainerv1.h>
#include <calobase/TowerInfoContainerv2.h>
#include <calobase/TowerInfoContainerv3.h>
#include <calobase/TowerInfoContainerv4.h>
#include <calobase/TowerInfoDefs.h>
#include <globalvertex/MbdVertex.h>
#include <globalvertex/MbdVertexv1.h>
using namespace std;
static const float radius_EM = 93.5;
static const float radius_OH = 225.87;

//____________________________________________________________________________..
zfinder::zfinder(const std::string &name, const int debug, const bool usez, const bool setz):
  SubsysReco(name)
{
  _setz = setz;
  _usez = usez;
  _name = name;
  _debug = debug;
  _nprocessed = 0;
}

//____________________________________________________________________________..
zfinder::~zfinder()
{

}

//____________________________________________________________________________..
int zfinder::Init(PHCompositeNode *topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int zfinder::InitRun(PHCompositeNode *topNode)
{
  if(_debug > 1) cout << "Initializing!" << endl;
  return Fun4AllReturnCodes::EVENT_OK;
}


float zfinder::new_eta(int channel, TowerInfoContainer* towers, RawTowerGeomContainer* geom, RawTowerDefs::CalorimeterId caloID, float testz)
{
  int key = towers->encode_key(channel);
  const RawTowerDefs::keytype geomkey = RawTowerDefs::encode_towerid(caloID, towers->getTowerEtaBin(key), towers->getTowerPhiBin(key));

  RawTowerGeom* tower_geom = geom->get_tower_geometry(geomkey);
  float oldeta = tower_geom->get_eta();
  
  float radius = (caloID==RawTowerDefs::CalorimeterId::HCALIN?radius_EM:radius_OH);
  float towerz = radius/(tan(2*atan(exp(oldeta))));
  float newz = towerz + testz;
  float newTheta = atan2(radius,newz);
  float neweta = -log(tan(0.5*newTheta));
  
  return neweta;
}

float get_dphi(float phi1, float phi2)
{
  float dphi = abs(phi1-phi2);
  if(dphi > M_PI) dphi = 2*M_PI - dphi;
  return dphi;
}

int zfinder::process_event(PHCompositeNode *topNode)
{

  

  if(_debug > 1) cout << endl << endl << endl << "zfinder: Beginning event processing" << endl;
  if(_nprocessed % 1000 == 0) cout << "processing event " << _nprocessed << endl;

  ++_nprocessed;
  GlobalVertexMap *globalmap = findNode::getClass<GlobalVertexMap>(topNode, "GlobalVertexMap");
  _zvtx = -300;
  JetContainer *jetcon = findNode::getClass<JetContainerv1>(topNode, "zzjets06");
  TowerInfoContainer *towers[3];
  towers[0] = findNode::getClass<TowerInfoContainer>(topNode, (!_usez)?"TOWERINFO_CALIB_CEMC_RETOWER":"TOWERINFO_CALIB_CEMC");
  towers[1] = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_HCALIN");
  towers[2] = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_HCALOUT");

  RawTowerGeomContainer *geom[3];
  geom[0] = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_CEMC");
  geom[1] = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALIN");
  geom[2] = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALOUT");


  float maxE[2] = {0};
  float oppmaxE[2] = {0};
  float maxz[2];
  float oppmaxz[2];
  float maxphi[2];
  float oppmaxphi[2];
  if(_usez)
    {
      if(_debug > 4) cout << "using z" << endl;
      for(int i=0; i<3; i+=2)
	{
	  RawTowerDefs::CalorimeterId caloID = (i==0)?RawTowerDefs::CalorimeterId::CEMC:RawTowerDefs::CalorimeterId::HCALOUT;
	  for(int j=0; j<((i==0)?24576:1536); ++j)
	    {
	      TowerInfo* tower = towers[i]->get_tower_at_channel(j);
	      if(tower->get_energy() > maxE[i/2])
		{
		  maxE[i/2] = tower->get_energy();
		  int key = towers[i]->encode_key(j);
		  const RawTowerDefs::keytype geomkey = RawTowerDefs::encode_towerid(caloID, towers[i]->getTowerEtaBin(key), towers[i]->getTowerPhiBin(key));

		  RawTowerGeom* tower_geom = geom[i]->get_tower_geometry(geomkey);
		  maxz[i/2] = tower_geom->get_center_z();
		  maxphi[i/2] = tower_geom->get_phi();
		}
	    }
	}
      for(int i=0; i<3; i+=2)
	{
	  RawTowerDefs::CalorimeterId caloID = (i==0)?RawTowerDefs::CalorimeterId::CEMC:RawTowerDefs::CalorimeterId::HCALOUT;
	  for(int j=0; j<((i==0)?24576:1536); ++j)
	    {
	      int key = towers[i]->encode_key(j);
	      const RawTowerDefs::keytype geomkey = RawTowerDefs::encode_towerid(caloID, towers[i]->getTowerEtaBin(key), towers[i]->getTowerPhiBin(key));
	      RawTowerGeom* tower_geom = geom[i]->get_tower_geometry(geomkey);
	      if(get_dphi(tower_geom->get_phi(), maxphi[i/2]) > 3*M_PI/4)
		{
		  TowerInfo* tower = towers[i]->get_tower_at_channel(j);
		  if(tower->get_energy() > maxE[i/2])
		    {
		      oppmaxE[i/2] = tower->get_energy();
		      oppmaxz[i/2] = tower_geom->get_center_z();
		      oppmaxphi[i/2] = tower_geom->get_phi();
		    }
		}
	    }
	}

      if(get_dphi(maxphi[0],oppmaxphi[0]) > M_PI/2)
	{
	  float tempE = oppmaxE[0];
	  float tempz = oppmaxz[0];
	  float tempphi = oppmaxphi[0];
	  oppmaxE[0] = oppmaxE[1];
	  oppmaxz[0] = oppmaxz[1];
	  oppmaxphi[0] = oppmaxphi[1];
	  oppmaxE[1] = tempE;
	  oppmaxz[1] = tempz;
	  oppmaxphi[1] = tempphi;
	}

      float avgz[2] = {0};
      float oppavgz[2] = {0};
      float sumE[2] = {0};
      float oppsumE[2] = {0};
      for(int i=0; i<3; i+=2)
	{
	  RawTowerDefs::CalorimeterId caloID = (i==0)?RawTowerDefs::CalorimeterId::CEMC:RawTowerDefs::CalorimeterId::HCALOUT;
	  for(int j=0; j<((i==0)?24576:1536); ++j)
	    {
	      int key = towers[i]->encode_key(j);
	      const RawTowerDefs::keytype geomkey = RawTowerDefs::encode_towerid(caloID, towers[i]->getTowerEtaBin(key), towers[i]->getTowerPhiBin(key));
	      RawTowerGeom* tower_geom = geom[i]->get_tower_geometry(geomkey);
	      if(get_dphi(tower_geom->get_phi(), maxphi[i/2]) < M_PI/4 && abs(maxz[i/2] - tower_geom->get_center_z()) < 100)
		{
		  TowerInfo* tower = towers[i]->get_tower_at_channel(j);
		  avgz[i/2] += tower->get_energy()*tower_geom->get_center_z();
		  sumE[i/2] += tower->get_energy();
		}
	      if(get_dphi(tower_geom->get_phi(), oppmaxphi[i/2]) < M_PI/4 && abs(oppmaxz[i/2] - tower_geom->get_center_z()) < 100)
		{
		  TowerInfo* tower = towers[i]->get_tower_at_channel(j);
		  oppavgz[i/2] += tower->get_energy()*tower_geom->get_center_z();
		  oppsumE[i/2] += tower->get_energy();
		}
	    }
	}

      if(_debug > 3) cout << avgz[0] << " " << avgz[1] << " " << sumE[0] << " " << sumE[1] << endl;

      for(int i=0; i<2; ++i)
	{
	  avgz[i] /= sumE[i];
	  oppavgz[i] /= oppsumE[i];
	}

      float slope[2] = {(radius_EM-radius_OH)/(avgz[0]-avgz[1]),(radius_EM-radius_OH)/(oppavgz[0]-oppavgz[1])};

      float zvals[2] = {avgz[0]-radius_EM/slope[0],oppavgz[0]-radius_EM/slope[1]};

      _zvtx = (zvals[0]+zvals[1])/2;
      
    }
  
  const int nz = 61;
  const int njet = 2;
  Jet* jets[njet];
  float jpt[njet] = {0};
  float jemsum[njet] = {0};
  float johsum[njet] = {0};
  float jemeta[njet] = {0};
  float joheta[njet] = {0};
  if(jetcon && !_usez)
    {
      int tocheck = jetcon->size();
      if(_debug > 2) cout << "Found " << tocheck << " jets to check..." << endl;
      for(int i=0; i<tocheck; ++i)
        {
          Jet *jet = jetcon->get_jet(i);
          if(jet)
            {
	      float pt = jet->get_pt();
	      if(pt < 10) continue;
	      if(pt > jpt[0])
		{
		  jets[0] = jet;
		  jpt[0] = pt;
		  jpt[1] = jpt[0];
		  jets[1] = jets[0];
		}
	      else if(pt > jpt[1])
		{
		  jets[1] = jet;
		  jpt[1] = pt;
		}
	    }
	  else
	    {
	      continue;
	    }
	}      
    }
  else
    {
      if(_debug > 0) cout << "no jets" << endl;
      if(!jetcon)
	{
	  return Fun4AllReturnCodes::ABORTEVENT;
	}
    }

  if(jpt[0] == 0 && _debug > 1 && !_usez)
    {
      cout << "NO JETS > 10 GeV!" << endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
	

  float metric = 999999999;
  if(!_usez)
    {
      for(int i=0; i<nz; ++i)
	{
	  float testz = -300+i*5;
	  float testmetric = 0;
	  for(int j=0; j<njet; ++j)
	    {
	      if(jpt[j] == 0) continue;
	      jemsum[j] = 0;
	      johsum[j] = 0;
	      jemeta[j] = 0;
	      joheta[j] = 0;
	      for(auto comp: jets[j]->get_comp_vec())
		{
		  if(comp.first==5 || comp.first == 26) continue;
		  unsigned int channel = comp.second;
		  if(comp.first==7 || comp.first == 27)
		    {
		      TowerInfo* tower = towers[2]->get_tower_at_channel(channel);
		      if(tower->get_energy() < 0.1) continue;
		      johsum[j] += tower->get_energy();
		      float neweta = new_eta(channel, towers[2], geom[2], RawTowerDefs::CalorimeterId::HCALOUT, testz);
		      joheta[j] += neweta*tower->get_energy();
		    }
		  if(comp.first == 13 || comp.first == 28 || comp.first == 25)
		    {
		      TowerInfo* tower = towers[0]->get_tower_at_channel(channel);
		      if(tower->get_energy() < 0.1) continue;
		      jemsum[j] += tower->get_energy();
		      float neweta = new_eta(channel, towers[0], geom[1], RawTowerDefs::CalorimeterId::HCALIN, testz);
		      jemeta[j] += neweta*tower->get_energy();
		    }
		}
	      if(_debug > 3) cout << jemeta[j] << " " << jemsum[j] << " : " << joheta[j] << " " << johsum[j] << endl;
	      jemeta[j] /= jemsum[j];
	      joheta[j] /= johsum[j];
	      testmetric += pow(jemeta[j]-joheta[j],2);
	    }
	  if(_debug > 3) cout << "metric: " << testmetric << endl;
	  if(testmetric < metric)
	    {
	      metric = testmetric;
	      _zvtx = testz;
	    }
	}
    }

  if(_debug > 2) cout << "optimal z: " << _zvtx << endl;
  if(_debug > 2) cout << globalmap->size() << endl;
  if(_setz)
    {
      MbdVertex *vertex = new MbdVertexv1();
      GlobalVertex* gvtx = globalmap->get(0);
      vertex->set_z(_zvtx);
      if(gvtx) gvtx->clone_insert_vtx(GlobalVertex::VTXTYPE::CALO,vertex);
      else
	{
	  gvtx = new GlobalVertexv2();
	  gvtx->clone_insert_vtx(GlobalVertex::VTXTYPE::CALO,vertex);
	  globalmap->insert(gvtx);
	}
      if(_debug > 2) cout << globalmap->size() << endl;
    }
  if(_debug > 3) cout << "end event" << endl;
  return Fun4AllReturnCodes::EVENT_OK;
    
}
//____________________________________________________________________________..
int zfinder::ResetEvent(PHCompositeNode *topNode)
{
  if (Verbosity() > 0)
    {
      std::cout << "zfinder::ResetEvent(PHCompositeNode *topNode) Resetting internal structures, prepare for next event" << std::endl;
    }
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int zfinder::EndRun(const int runnumber)
{
  if (Verbosity() > 0)
    {
      std::cout << "zfinder::EndRun(const int runnumber) Ending Run for Run " << runnumber << std::endl;
    }
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int zfinder::End(PHCompositeNode *topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int zfinder::Reset(PHCompositeNode *topNode)
{
  if (Verbosity() > 0)
    {
      std::cout << "zfinder::Reset(PHCompositeNode *topNode) being Reset" << std::endl;
    }
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
void zfinder::Print(const std::string &what) const
{
  std::cout << "zfinder::Print(const std::string &what) const Printing info for " << what << std::endl;
}
