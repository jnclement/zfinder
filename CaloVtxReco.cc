#include "CaloVtxReco.h"
#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>
#include <jetbase/JetContainerv1.h>
#include <jetbase/Jet.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/RawTowerGeomContainer_Cylinderv1.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfoContainerv1.h>
#include <calobase/TowerInfoContainerv2.h>
#include <calobase/TowerInfoContainerv3.h>
#include <calobase/TowerInfoContainerv4.h>
#include <calobase/TowerInfoDefs.h>
#include <globalvertex/CaloVertexv1.h>
#include <globalvertex/CaloVertexMapv1.h>
using namespace std;
static const float radius_EM = 93.5;
static const float radius_OH = 225.87;

//____________________________________________________________________________..
CaloVtxReco::CaloVtxReco(const std::string &name, const int debug, const bool usez, const bool setz):
  SubsysReco(name)
{
  _setz = setz;
  _usez = usez;
  _name = name;
  _debug = debug;
  _nprocessed = 0;
}

//____________________________________________________________________________..
CaloVtxReco::~CaloVtxReco()
{

}

int CaloVtxReco::createNodes(PHCompositeNode *topNode)
{
  PHNodeIterator iter(topNode);
  PHCompositeNode *dstNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "DST"));
  if (!dstNode)
  {
    std::cout << PHWHERE << "DST Node missing doing nothing" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  PHCompositeNode *runNode = dynamic_cast<PHCompositeNode *>(iter.findFirst("PHCompositeNode", "RUN"));
  if (!runNode)
  {
    std::cout << PHWHERE << "RUN Node missing doing nothing" << std::endl;
    return Fun4AllReturnCodes::ABORTRUN;
  }

  PHNodeIterator dstiter(dstNode);
  
  PHCompositeNode *globalNode = dynamic_cast<PHCompositeNode *>(dstiter.findFirst("PHCompositeNode", "GLOBAL"));
  if (!globalNode)
  {
    globalNode = new PHCompositeNode("GLOBAL");
    dstNode->addNode(globalNode);
  }

  _calovtxmap = findNode::getClass<CaloVertexMap>(globalNode, "CaloVertexMap");
  if (!_calovtxmap)
  {
    _calovtxmap = new CaloVertexMapv1();
    PHIODataNode<PHObject> *VertexMapNode = new PHIODataNode<PHObject>(_calovtxmap, "CaloVertexMap", "PHObject");
    globalNode->addNode(VertexMapNode);
  }

  return Fun4AllReturnCodes::EVENT_OK;
}


//____________________________________________________________________________..
int CaloVtxReco::Init(PHCompositeNode *topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int CaloVtxReco::InitRun(PHCompositeNode *topNode)
{
  if(_debug > 1) cout << "Initializing!" << endl;
  if(createNodes(topNode) == Fun4AllReturnCodes::ABORTEVENT)
    {
      return Fun4AllReturnCodes::ABORTEVENT;
    }
  return Fun4AllReturnCodes::EVENT_OK;
}


float CaloVtxReco::new_eta(int channel, TowerInfoContainer* towers, RawTowerGeomContainer* geom, RawTowerDefs::CalorimeterId caloID, float testz)
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

int CaloVtxReco::process_event(PHCompositeNode *topNode)
{

  

  if(_debug > 1) cout << endl << endl << endl << "CaloVtxReco: Beginning event processing" << endl;
  if(_nprocessed % 1000 == 0) cout << "processing event " << _nprocessed << endl;

  ++_nprocessed;
  _zvtx = -9999;
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
  const int nz = 121;
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
	      if(pt < 5) continue;
	      if(pt > jpt[0])
		{
		  jpt[1] = jpt[0];
		  jets[1] = jets[0];
		  jets[0] = jet;
		  jpt[0] = pt;
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
      if(!jetcon && !_usez)
	{
	  return Fun4AllReturnCodes::ABORTEVENT;
	}
    }

  if(jpt[0] == 0 && !_usez)
    {
      if(_debug > 2) cout << "NO JETS > 5 GeV!" << endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
	

  float metric = FLT_MAX;
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
		  if(_debug > 6) cout << comp.first << " ";
		  if(comp.first==5 || comp.first == 26) continue;
		  unsigned int channel = comp.second;
		  if(comp.first==7 || comp.first == 27)
		    {
		      TowerInfo* tower = towers[2]->get_tower_at_channel(channel);
		      if(_debug > 6) cout << "towerE: " << tower->get_energy() << " " << endl;
		      if(tower->get_energy() < 0.1) continue;
		      if(_debug > 6) cout << "good tower ";
		      johsum[j] += tower->get_energy();
		      float neweta = new_eta(channel, towers[2], geom[2], RawTowerDefs::CalorimeterId::HCALOUT, testz);
		      joheta[j] += neweta*tower->get_energy();
		    }
		  if(comp.first == 13 || comp.first == 28 || comp.first == 25)
		    {
		      TowerInfo* tower = towers[0]->get_tower_at_channel(channel);
		      if(_debug > 6) cout << "towerE: " << tower->get_energy() << " " << endl;
		      if(tower->get_energy() < 0.1) continue;
		      if(_debug > 6) cout << "good tower ";
		      jemsum[j] += tower->get_energy();
		      float neweta = new_eta(channel, towers[0], geom[1], RawTowerDefs::CalorimeterId::HCALIN, testz);
		      jemeta[j] += neweta*tower->get_energy();
		    }
		}
	      if(_debug > 6) cout << endl;
	      if(_debug > 3) cout << jemeta[j] << " " << jemsum[j] << " : " << joheta[j] << " " << johsum[j] << endl;
	      jemeta[j] /= jemsum[j];
	      joheta[j] /= johsum[j];
	      if((jemsum[j] == 0 || johsum[j] == 0) && _debug > 1) cout << "zero E sum in at least one calo for a jet" << endl;
	      testmetric += pow(jemeta[j]-joheta[j],2);
	    }
	  if(_debug > 3) cout << "metric: " << testmetric << endl;
	  if(testmetric < metric && testmetric != 0)
	    {
	      metric = testmetric;
	      _zvtx = testz;
	    }
	}
    }

  if(_debug > 2) cout << "optimal z: " << _zvtx << endl;
  if(_setz)
    {
      CaloVertex *vertex = new CaloVertexv1();
      if(!_usez) _zvtx *= 1.406; //calibration factor from simulation
      vertex->set_z(_zvtx);
      _calovtxmap->insert(vertex);
    }
  if(_debug > 3) cout << "end event" << endl;
  return Fun4AllReturnCodes::EVENT_OK;
    
}
//____________________________________________________________________________..
int CaloVtxReco::ResetEvent(PHCompositeNode *topNode)
{
  if (Verbosity() > 0)
    {
      std::cout << "CaloVtxReco::ResetEvent(PHCompositeNode *topNode) Resetting internal structures, prepare for next event" << std::endl;
    }
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int CaloVtxReco::EndRun(const int runnumber)
{
  if (Verbosity() > 0)
    {
      std::cout << "CaloVtxReco::EndRun(const int runnumber) Ending Run for Run " << runnumber << std::endl;
    }
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int CaloVtxReco::End(PHCompositeNode *topNode)
{
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int CaloVtxReco::Reset(PHCompositeNode *topNode)
{
  if (Verbosity() > 0)
    {
      std::cout << "CaloVtxReco::Reset(PHCompositeNode *topNode) being Reset" << std::endl;
    }
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
void CaloVtxReco::Print(const std::string &what) const
{
  std::cout << "CaloVtxReco::Print(const std::string &what) const Printing info for " << what << std::endl;
}
