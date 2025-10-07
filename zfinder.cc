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
zfinder::zfinder(const std::string &name, const int debug):
  SubsysReco(name)
{
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


int zfinder::process_event(PHCompositeNode *topNode)
{

  

  if(_debug > 1) cout << endl << endl << endl << "zfinder: Beginning event processing" << endl;
  if(_nprocessed % 1000 == 0) cout << "processing event " << _nprocessed << endl;

  ++_nprocessed;
  GlobalVertexMap *globalmap = findNode::getClass<GlobalVertexMap>(topNode, "GlobalVertexMap");
  _zvtx = -150;
  JetContainer *jetcon = findNode::getClass<JetContainerv1>(topNode, "zzjets");
  TowerInfoContainer *towers[3];
  towers[0] = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_CEMC_RETOWER");
  towers[1] = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_HCALIN");
  towers[2] = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_HCALOUT");

  RawTowerGeomContainer *geom[3];
  geom[0] = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_CEMC");
  geom[1] = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALIN");
  geom[2] = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALOUT");
  
  const int nz = 61;
  const int njet = 2;
  Jet* jets[njet];
  float jpt[njet] = {0};
  float jemsum[njet] = {0};
  float johsum[njet] = {0};
  float jemeta[njet] = {0};
  float joheta[njet] = {0};
  if(jetcon)
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
      return Fun4AllReturnCodes::ABORTEVENT;
    }

  if(jpt[0] == 0 && _debug > 1)
    {
      cout << "NO JETS > 10 GeV!" << endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
	

  float metric = 999999999;
  
  for(int i=0; i<nz; ++i)
    {
      float testz = -150+i*5;
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
		  johsum[j] += tower->get_energy();
		  float neweta = new_eta(channel, towers[2], geom[2], RawTowerDefs::CalorimeterId::HCALOUT, testz);
		  joheta[j] += neweta*tower->get_energy();
		}
	      if(comp.first == 13 || comp.first == 28 || comp.first == 25)
		{
		  TowerInfo* tower = towers[0]->get_tower_at_channel(channel);
		  jemsum[j] += tower->get_energy();
		  float neweta = new_eta(channel, towers[0], geom[1], RawTowerDefs::CalorimeterId::HCALIN, testz);
		  jemeta[j] += neweta*tower->get_energy();
		}
	    }
	  jemeta[j] /= jemsum[j];
	  joheta[j] /= johsum[j];
	  testmetric += pow(jemeta[j]-joheta[j],2);
	}
      if(testmetric < metric)
	{
	  metric = testmetric;
	  _zvtx = testz;
	}
    }

  if(_debug > 2) cout << "optimal z: " << _zvtx << endl;
  if(_debug > 2) cout << globalmap->size() << endl;
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
