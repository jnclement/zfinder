#include <iomanip>
#include "zfinder.h"
#include <calobase/RawTower.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerGeom.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/RawTowerGeomContainer_Cylinderv1.h>
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfoContainerv1.h>
#include <calobase/TowerInfoContainerv2.h>
#include <calobase/TowerInfoContainerv3.h>
#include <calobase/TowerInfoContainerv4.h>
#include <globalvertex/GlobalVertexMap.h>
#include <globalvertex/GlobalVertexMapv1.h>
#include <globalvertex/GlobalVertex.h>
#include <mbd/MbdPmtContainer.h>
#include <mbd/MbdPmtContainerV1.h>
#include <fun4all/Fun4AllReturnCodes.h>
#include <phool/PHCompositeNode.h>
#include <phool/PHRandomSeed.h>
#include <phool/getClass.h>
#include <globalvertex/MbdVertexMap.h>
#include <globalvertex/MbdVertexMapv1.h>
#include <globalvertex/MbdVertex.h>
#include <g4main/PHG4TruthInfoContainer.h>
#include <g4main/PHG4Particle.h>
#include <mbd/MbdPmtHit.h>
#include <jetbackground/TowerBackgroundv1.h>
#include <cmath>
#include <mbd/MbdOut.h>
#include <jetbase/FastJetAlgo.h>
#include <jetbase/JetReco.h>
#include <jetbase/TowerJetInput.h>
#include <g4jets/TruthJetInput.h>
#include <jetbase/JetContainerv1.h>
#include <jetbase/Jet.h>
#include <calobase/RawTowerv1.h>
#include <jetbackground/CopyAndSubtractJets.h>
#include <jetbackground/DetermineTowerBackground.h>
#include <jetbackground/FastJetAlgoSub.h>
#include <jetbackground/RetowerCEMC.h>
#include <jetbackground/SubtractTowers.h>
#include <jetbackground/SubtractTowersCS.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_rng.h>  // for gsl_rng_uniform_pos
#include <phool/recoConsts.h>
#include <iostream>
#include <centrality/CentralityInfo.h>
#include <calotrigger/MinimumBiasInfo.h>
#include <TH2D.h>
#include <TColor.h>
#include <TMarker.h>
#include <TPad.h>
#include <TCanvas.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TLine.h>
#include <TH1D.h>
#include <calotrigger/MinimumBiasInfo.h>
#include <calotrigger/MinimumBiasInfov1.h>
#include <calotrigger/MinimumBiasClassifier.h>
#include <ffarawobjects/Gl1Packetv3.h>
#include <TLorentzVector.h>
#include <ffaobjects/EventHeader.h>
#include <pdbcalbase/PdbParameterMap.h>
#include <calobase/TowerInfoDefs.h>
#include <ffarawobjects/CaloPacketv1.h>
using namespace std;
static const float radius_EM = 93.5;
static const float minz_EM = -130.23;
static const float maxz_EM = 130.23;

static const float radius_IH = 127.503;
static const float minz_IH = -170.299;
static const float maxz_IH = 170.299;

static const float radius_OH = 225.87;
static const float minz_OH = -301.683;
static const float maxz_OH = 301.683;
float get_emcal_mineta_zcorrected(float zvertex) {
  float z = minz_EM - zvertex;
  float eta_zcorrected = asinh(z / (float)radius_EM);
  return eta_zcorrected;
}

float get_emcal_maxeta_zcorrected(float zvertex) {
  float z = maxz_EM - zvertex;
  float eta_zcorrected = asinh(z / (float)radius_EM);
  return eta_zcorrected;
}

float get_ihcal_mineta_zcorrected(float zvertex) {
  float z = minz_IH - zvertex;
  float eta_zcorrected = asinh(z / (float)radius_IH);
  return eta_zcorrected;
}

float get_ihcal_maxeta_zcorrected(float zvertex) {
  float z = maxz_IH - zvertex;
  float eta_zcorrected = asinh(z / (float)radius_IH);
  return eta_zcorrected;
}

float get_ohcal_mineta_zcorrected(float zvertex) {
  float z = minz_OH - zvertex;
  float eta_zcorrected = asinh(z / (float)radius_OH);
  return eta_zcorrected;
}

float get_ohcal_maxeta_zcorrected(float zvertex) {
  float z = maxz_OH - zvertex;
  float eta_zcorrected = asinh(z / (float)radius_OH);
  return eta_zcorrected;
}

bool check_bad_jet_eta(float jet_eta, float zertex, float jet_radius) {
  float emcal_mineta = get_emcal_mineta_zcorrected(zertex);
  float emcal_maxeta = get_emcal_maxeta_zcorrected(zertex);
  float ihcal_mineta = get_ihcal_mineta_zcorrected(zertex);
  float ihcal_maxeta = get_ihcal_maxeta_zcorrected(zertex);
  float ohcal_mineta = get_ohcal_mineta_zcorrected(zertex);
  float ohcal_maxeta = get_ohcal_maxeta_zcorrected(zertex);
  float minlimit = emcal_mineta;
  if (ihcal_mineta > minlimit) minlimit = ihcal_mineta;
  if (ohcal_mineta > minlimit) minlimit = ohcal_mineta;
  float maxlimit = emcal_maxeta;
  if (ihcal_maxeta < maxlimit) maxlimit = ihcal_maxeta;
  if (ohcal_maxeta < maxlimit) maxlimit = ohcal_maxeta;
  minlimit += jet_radius;
  maxlimit -= jet_radius;
  return jet_eta < minlimit || jet_eta > maxlimit;
}

void drawText(const char *text, float xp, float yp, bool isRightAlign=0, int textColor=kBlack, double textSize=0.04, int textFont = 42, bool isNDC=true){
  // when textfont 42, textSize=0.04                                                                                                                         
  // when textfont 43, textSize=18                                                                                                                           
  TLatex *tex = new TLatex(xp,yp,text);
  tex->SetTextFont(textFont);
  //   if(bold)tex->SetTextFont(43);                                                                                                                         
  tex->SetTextSize(textSize);
  tex->SetTextColor(textColor);
  tex->SetLineWidth(1);
  if(isNDC) tex->SetNDC();
  if(isRightAlign) tex->SetTextAlign(31);
  tex->Draw();
  //delete tex;
}

void sqrt_snn_text(float xp = 0.7, float yp = 0.8, bool isRightAlign=0, double textsize = 0.04)
{
  drawText("#sqrt{S_{NN}} = 200 GeV",xp,yp,isRightAlign,kBlack,textsize);
}

void sphenixtext(float xpos = 0.7, float ypos = 0.96, int ra = 0, float textsize = 0.04)
{
  drawText("#bf{#it{sPHENIX}} Internal", xpos, ypos, ra, kBlack, textsize);
}


float zfinder::bintoeta_hc(int etabin)
{
  return (2.2*etabin)/24 - 1.1;
}
float zfinder::bintophi_hc(int phibin)
{
  return (2*M_PI*phibin)/64;
}

float bintoeta_em(int etabin)
{
  return (2.2*etabin)/96-1.1;
}

float bintophi_em(int phibin)
{
  return (2*M_PI*phibin)/256;
}

void zfinder::drawCalo(TowerInfoContainer** towers, float* jet_e, float* jet_et, float* jet_ph, int jet_n, float jet_ecc, float jet_lfrac, RawTowerGeomContainer** geom, float zvtx, int failscut, int runnum, int evtnum, float frcoh, float frcem, float emtot, float ohtot, float maxJetE)
{
  std::stringstream full_stream;
  full_stream << std::fixed << std::setprecision(3) << "Run " << runnum << ", event " << evtnum << ", EM fraction: " << frcem << ", OH fraction: " << frcoh << ", z_{vtx} = " <<  std::setprecision(1) << zvtx << " ";
  std::string full_string = full_stream.str();

  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);
  gStyle->SetNumberContours(50);
  int ncircle = 64;
  const int ncol = 50;
  const int nstp = 3;
  double red[nstp] = {0.,1.,1.};// = {1, 1, 1, 1, 1, 1, 1, 1};
  double grn[nstp] = {0.,1.,0.};// = {1, 0.875, 0.75, 0.625, 0.5, 0.375, 0.25, 0.125, 0};
  double blu[nstp] = {1.,1.,0.};// = {1, 0.875, 0.75, 0.625, 0.5, 0.375, 0.25, 0.125, 0};
  double stp[nstp] = {0,2./27,1};;// = {0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875, 1.0};
  

  
  TCanvas* c = new TCanvas("","",1900,1000);
  c->Divide(4,1,0,0.1);
  
  TH2D* event_sum = new TH2D("event_sum","Calorimeter Sum",24,-0.5,23.5,64,-0.5,63.5);
  TH2D* event_disrt[3];
  for(int i=0; i<3; ++i)
    {
      int nbinx = (i==0?96:24);
      int nbiny = (i==0?256:64);
      event_disrt[i] = new TH2D(("event_display_rt"+to_string(i)).c_str(),"",nbinx,-0.5,nbinx-0.5,nbiny,-0.5,nbiny-0.5);

    }
  TColor::CreateGradientColorTable(nstp, stp, red, grn, blu, ncol);
  event_disrt[0]->GetXaxis()->SetTitle("EMCal #eta");
  event_disrt[0]->GetYaxis()->SetTitle("EMCal #phi");
  event_disrt[1]->GetXaxis()->SetTitle("IHCal #eta");
  event_disrt[1]->GetYaxis()->SetTitle("IHCal #phi");
  event_disrt[2]->GetXaxis()->SetTitle("OHCal #eta");
  event_disrt[2]->GetYaxis()->SetTitle("OHCal #phi");
  for(int i=0; i<3; ++i)
    {
      event_disrt[i]->GetZaxis()->SetTitleOffset(1);
      event_disrt[i]->GetYaxis()->SetTitleOffset(1);
      event_disrt[i]->GetZaxis()->SetTitle("Uncorrected Tower Energy [GeV]");
      event_disrt[i]->GetZaxis()->SetRangeUser(-2,25);
      event_disrt[i]->GetXaxis()->SetNdivisions(4,kFALSE);
      event_disrt[i]->GetXaxis()->SetTitleSize(0.04);
      event_disrt[i]->GetYaxis()->SetTitleSize(0.04);
      event_disrt[i]->GetZaxis()->SetTitleSize(0.04);
      event_disrt[i]->GetXaxis()->SetTitleOffset(1);
      event_disrt[i]->GetXaxis()->SetLabelSize(0.04);
      event_disrt[i]->GetYaxis()->SetLabelSize(0.04);
      event_disrt[i]->GetZaxis()->SetLabelSize(0.04);
      event_disrt[i]->GetXaxis()->SetLabelOffset(0.02);
    }
  event_sum->GetXaxis()->SetTitle("Calo Sum #eta");
  event_sum->GetYaxis()->SetTitle("Calo Sum #phi");
  event_sum->GetZaxis()->SetTitleOffset(1);
  event_sum->GetYaxis()->SetTitleOffset(1);
  event_sum->GetZaxis()->SetTitle("Uncorrected Tower Energy [GeV]");
  event_sum->GetZaxis()->SetRangeUser(-2,25);
  event_sum->GetXaxis()->SetNdivisions(4,kFALSE);
  event_sum->GetXaxis()->SetTitleSize(0.04);
  event_sum->GetYaxis()->SetTitleSize(0.04);
  event_sum->GetZaxis()->SetTitleSize(0.04);
  event_sum->GetXaxis()->SetTitleOffset(1);
  event_sum->GetXaxis()->SetLabelSize(0.04);
  event_sum->GetYaxis()->SetLabelSize(0.04);
  event_sum->GetZaxis()->SetLabelSize(0.04);
  event_sum->GetXaxis()->SetLabelOffset(0.02);
  event_sum->Reset();
  for(int j=0; j<3; ++j)
    {
      event_disrt[j]->Reset();
      for(int k=0; k<(j==0?24576:1536); ++k)
	{
	  TowerInfo* tower = towers[j]->get_tower_at_channel(k);
	  int key = towers[j]->encode_key(k);
	  float eta = towers[j]->getTowerEtaBin(key);
	  float phi = towers[j]->getTowerPhiBin(key);
	  event_disrt[j]->Fill(eta,phi,tower->get_energy());
	  if(j==0)event_sum->Fill(eta/4,phi/4,tower->get_energy());
	  else event_sum->Fill(eta,phi,tower->get_energy());
	  /*
	  if(j==1)
	    {
	      
	      const RawTowerDefs::keytype geomkey = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::HCALIN, towers[j]->getTowerEtaBin(key), towers[j]->getTowerPhiBin(key));
	      RawTowerGeom *tower_geom = geom[j]->get_tower_geometry(geomkey); //encode tower geometry
	      float newx = tower_geom->get_center_x();
	      float newy = tower_geom->get_center_y();
	      float newz = tower_geom->get_center_z() - zvtx;
	
	      eta = asinh(newz/sqrt(newx*newx+newy*newy));//getEtaFromBinEM()+0.012;
			  
	      phi = tower_geom->get_phi()+M_PI;//bintophi_hc(towers[j]->getTowerPhiBin(key))+0.048;
	    }
	  else if(j==2)
	    {
	      const RawTowerDefs::keytype geomkey = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::HCALOUT, towers[j]->getTowerEtaBin(key), towers[j]->getTowerPhiBin(key));
	      RawTowerGeom *tower_geom = geom[j]->get_tower_geometry(geomkey); //encode tower geometry
	      float newx = tower_geom->get_center_x();
	      float newy = tower_geom->get_center_y();
	      float newz = tower_geom->get_center_z() - zvtx;
	      eta = asinh(newz/sqrt(newx*newx+newy*newy));//getEtaFromBinEM()+0.012;
	      phi = tower_geom->get_phi()+M_PI;//
	    }
	  else
	    {
	      const RawTowerDefs::keytype geomkey = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::HCALIN, towers[j]->getTowerEtaBin(key), towers[j]->getTowerPhiBin(key));
	      RawTowerGeom *tower_geom = geom[1]->get_tower_geometry(geomkey); //encode tower geometry

	      float radius = 93.5;
	      float ihEta = tower_geom->get_eta();
	      float emZ = radius/(tan(2*atan(exp(-ihEta))));
	      float newz = emZ - zvtx;
	      float newTheta = atan2(radius,newz);
	      float towerEta = -log(tan(0.5*newTheta));
	      eta = towerEta;
	      phi = tower_geom->get_phi()+M_PI;//bintophi_em(towers[j]->getTowerPhiBin(key))+0.012;
	    }
	  */

	}


      c->cd(j+1);
      gPad->SetLogz(0);
      gPad->SetRightMargin(0.2);
      gPad->SetLeftMargin(0.2);
      gPad->SetTopMargin(0.05);
      event_disrt[j]->SetContour(ncol);
      event_disrt[j]->Draw("COLZ");

    }
  c->cd(4);
  gPad->SetLogz(0);            
  gPad->SetTopMargin(0.05);                                                       
  gPad->SetRightMargin(0.2);
  gPad->SetLeftMargin(0.2);
  event_sum->SetContour(ncol);
  event_sum->Draw("COLZ");
  c->cd(0);
  sphenixtext(0.96,0.96,1,0.04);

  std::stringstream ecc_stream;
  ecc_stream << std::fixed << std::setprecision(3) << jet_ecc;
  std::string ecc_string = ecc_stream.str();
  //drawText(("#epsilon_{max jet}="+ecc_string).c_str(),0.1,0.95);

  std::stringstream lfrac_stream;
  lfrac_stream << std::fixed << std::setprecision(3) << jet_lfrac;
  std::string lfrac_string = lfrac_stream.str();
  //drawText(("E_{max layer,jet}/E_{max jet}="+lfrac_string).c_str(),0.4,0.95);

  std::stringstream z_stream;
  z_stream << std::fixed << std::setprecision(1) << zvtx;
  std::string z_string = z_stream.str();
  //  drawText(("z_{vtx}="+z_string).c_str(),0.25,0.95);
  
  std::string fails = "";
  std::string whichcut = "";
  if(failscut == 0)
    {
      fails = "Passes dijet cut";
      whichcut = "dijet";
    }
  else if(failscut == 1)
    {
      fails = "Passes frac. cut";
      whichcut = "frac";
    }
  else
    {
      fails = "Passes both";
      whichcut = "both";
    }
  full_string += fails;
  
  //drawText(fails.c_str(),0.7,0.95);
  drawText(full_string.c_str(),0.05,0.95,0,kBlack,0.02);
  c->cd(4);
  maxJetE = 0;
  for(int k=0; k<jet_n; ++k)
    {
      if(maxJetE < jet_e[k])
	{
	  maxJetE = jet_e[k];
	}
      for(int l=0; l<ncircle; ++l)
	{
	  float eta = jet_et[k]+0.4*cos(2*l*M_PI/ncircle);
	  float phi = jet_ph[k]+0.4*sin(2*l*M_PI/ncircle);
	  phi += M_PI;
	  if(eta > 2.2 || eta < -2.2) continue;
	  if(phi > 2*M_PI) phi -= 2*M_PI;
	  if(phi < 0) phi += 2*M_PI;
	  TMarker* circlemarker = new TMarker(eta,phi,20);
	  circlemarker->SetMarkerSize(0.3);
	  circlemarker->SetMarkerColor(kBlue);
	  //circlemarker->Draw();
        }
      std::stringstream e_stream;
      e_stream << std::fixed << std::setprecision(2) << jet_e[k];
      std::string e_string = e_stream.str();
      drawText((e_string+" GeV").c_str(),12,((jet_ph[k]/*+(jet_ph[k]+M_PI>3.84?-0.53:0.43)*/+(jet_ph[k]<0?2*M_PI:0))/(2*M_PI))*64,/*(jet_et[k]>0?1:*/0/*)*/,kBlack,0.04,42,false);

    }
  //c->Update();
  string dirstring = "";
  for(int i=46; i<101; ++i)
    {
      if(maxJetE < i)
	{
	  dirstring = to_string(i-1)+"to"+to_string(i);
	  break;
	}
    }
  if(maxJetE > 100) dirstring = "gr100";
      
  c->SaveAs(("/sphenix/user/jocl/projects/run2024_earlydata/run/output/smg/candidate_"+to_string(runnum)+"_"+dirstring+"_"+_name+"_"+whichcut+"_"+to_string(cancount)+".png").c_str());
  cout << "Saved" << endl;
  cancount++;

  for(int i=0; i<3; ++i)
    {
      c->cd(i+1);
      gPad->SetLogz();
      event_disrt[i]->GetZaxis()->SetRangeUser(0.05,25);
      gPad->Update();
    }
  c->cd(4);
  gPad->SetLogz();
  event_sum->GetZaxis()->SetRangeUser(0.05,25);
  gPad->Update();
  c->SaveAs(("/sphenix/user/jocl/projects/run2024_earlydata/run/output/smg/candidate_"+to_string(runnum)+"_"+dirstring+"_"+_name+"_"+whichcut+"_"+to_string(cancount)+"_log.png").c_str());
  /*
  delete c;
  delete event_sum;
  for(int i=0; i<3; ++i)
    {
      delete event_disrt[i];
    }
  */
}

//____________________________________________________________________________..
zfinder::zfinder(const std::string &filename, const std::string &name, const int debug, const std::string &wfilename, const int dowf, const bool isdat, const int doall60, const int dotruthpar):
  SubsysReco(name), _cutParams(name)
{
  _dotruthpar = dotruthpar;
  _isdat = isdat;
  _name = name;
  _debug = debug;
  _filename = filename;
  _wfilename = wfilename;
  _nprocessed = 0;
  _dowf = dowf;
  jet_tree = new TTree("jet_tree","a persevering date tree");
  if(_dowf) _wft = new TTree("wft","a stupid waveform tree");
  _doall60 = doall60;
  if(_doall60 && _dowf) _wff = new TFile(_wfilename.c_str(),"RECREATE");
  if(_doall60 && _dowf) _wff->cd();
  if(_doall60 && _dowf) _wft->SetDirectory(_wff);
  //if(_doall60 || !_isdat)
  _f = new TFile(_filename.c_str(), "RECREATE");
  _f->cd();
  jet_tree->SetDirectory(_f);
}

//____________________________________________________________________________..
zfinder::~zfinder()
{

}

//____________________________________________________________________________..
int zfinder::Init(PHCompositeNode *topNode)
{
  
  if(_debug > 1) cout << "Begin init: " << endl;
  /*
  jet_tree->Branch("nBadChi2",&_nBadChi2,"nBadChi2/I");
  jet_tree->Branch("ecc",&_eccentricity,"ecc/F");
  jet_tree->Branch("frcoh",&_frcoh,"frcoh/F");
  jet_tree->Branch("frcem",&_frcem,"frcem/F");



  jet_tree->Branch("eta",&_eta,"eta/F");
  jet_tree->Branch("phi",&_phi,"phi/F");
  jet_tree->Branch("jet_ET",&_jet_ET,"jet_ET/F");
  jet_tree->Branch("dphi",&_dphi,"dphi/F");
  jet_tree->Branch("isdijet",&_isdijet,"isdijet/I");
  jet_tree->Branch("subjet_ET",&_subjet_ET,"subjet_ET/F");
  */
  jet_tree->Branch("zvtx",&_zvtx,"zvtx/F");
  //jet_tree->Branch("jetcompE",_jetcompE,"jetcompE[3][512]/F");
  //jet_tree->Branch("jetcompEta",_jetcompEta,"jetcompEta[3][512]/F");
  //jet_tree->Branch("jetcompPhi",_jetcompPhi,"jetcompPhi[3][512]/F");
  /*
  jet_tree->Branch("maxTowChi2",_maxTowChi2,"maxTowChi2[3]/F");
  jet_tree->Branch("maxTowE",&_maxTowE,"maxTowE/F");
  jet_tree->Branch("subTowE",&_subTowE,"subTowE/F");
  jet_tree->Branch("maxTowDiff",&_maxTowDiff,"maxTowDiff/F");
  jet_tree->Branch("maxETowChi2",&_maxETowChi2,"maxETowChi2/F");
  jet_tree->Branch("ohPhiBinMaxFrac",&_ohPhiBinMaxFrac,"ohPhiBinMaxFrac/F");
  jet_tree->Branch("maxETowIsZS",&_maxETowIsZS,"maxETowIsZS/I");
  jet_tree->Branch("maxETowChi2Det",&_maxETowChi2Det,"maxETowChi2Det/I");
  jet_tree->Branch("triggervec",&_triggervec,"triggervec/g");
  */
  jet_tree->Branch("bbfqavec",&_bbfqavec,"bbfqavec/i");
  /*
  jet_tree->Branch("elmbgvec",&_elmbgvec,"elmbgvec/i");
  */
  jet_tree->Branch("jet_n",&_jet_n,"jet_n/I");
  jet_tree->Branch("runnum",&_runnum,"runnum/I");
  jet_tree->Branch("evtnum",&_evtnum,"evtnum/I");
  jet_tree->Branch("failscut",&_failscut,"failscut/I");
  jet_tree->Branch("alljetfrcoh",_alljetfrcoh,"alljetfrcoh[jet_n]/F");
  jet_tree->Branch("alljetfrcem",_alljetfrcem,"alljetfrcem[jet_n]/F");
  jet_tree->Branch("jet_et",_jet_et,"jet_et[jet_n]/F");
  jet_tree->Branch("jet_pt",_jet_pt,"jet_pt[jet_n]/F");
  jet_tree->Branch("jet_t",_jet_t,"jet_t[jet_n]/F");
  //jet_tree->Branch("jet_t_em",_jet_t_em,"jet_t_em[jet_n]/F");
  //jet_tree->Branch("jet_t_ih",_jet_t_ih,"jet_t_ih[jet_n]/F");
  //jet_tree->Branch("jet_t_oh",_jet_t_oh,"jet_t_oh[jet_n]/F");
  jet_tree->Branch("jet_eta",_jet_eta,"jet_eta[jet_n]/F");
  jet_tree->Branch("jet_phi",_jet_phi,"jet_phi[jet_n]/F");

  if(!_isdat)
    {
  jet_tree->Branch("tjet_n",&_tjet_n,"tjet_n/I");
  jet_tree->Branch("tjet_e",_tjet_e,"tjet_e[tjet_n]/F");
  jet_tree->Branch("tjet_pt",_tjet_pt,"tjet_pt[tjet_n]/F");
  jet_tree->Branch("tjet_eta",_tjet_eta,"tjet_eta[tjet_n]/F");
  jet_tree->Branch("tjet_phi",_tjet_phi,"tjet_phi[tjet_n]/F");
    }
  /*
  jet_tree->Branch("emtow",_emtow,"emtow[96][256]/F");
  jet_tree->Branch("ihtow",_ihtow,"ihtow[24][64]/F");
  jet_tree->Branch("ohtow",_ohtow,"ohtow[24][64]/F");
  jet_tree->Branch("isbadem",_isbadem,"isbadem[96][256]/I");
  jet_tree->Branch("isbadih",_isbadih,"isbadih[24][64]/I");
  jet_tree->Branch("isbadoh",_isbadoh,"isbadoh[24][64]/I");
  jet_tree->Branch("ishotem",_ishotem,"ishotem[96][256]/I");
  jet_tree->Branch("ishotih",_ishotih,"ishotih[24][64]/I");
  jet_tree->Branch("ishotoh",_ishotoh,"ishotoh[24][64]/I");
  jet_tree->Branch("nocalem",_nocalem,"nocalem[96][256]/I");
  jet_tree->Branch("nocalih",_nocalih,"nocalih[24][64]/I");
  jet_tree->Branch("nocaloh",_nocaloh,"nocaloh[24][64]/I");
  jet_tree->Branch("jconem",_jconem,"jconem[24][64]/F");
  jet_tree->Branch("jconih",_jconih,"jconih[24][64]/F");
  jet_tree->Branch("jconoh",_jconoh,"jconoh[24][64]/F");
  jet_tree->Branch("isblt",&_isbadlive,"isblt/I");
  jet_tree->Branch("chi2em",_chi2em,"chi2em[96][256]/F");
  jet_tree->Branch("chi2ih",_chi2ih,"chi2ih[24][64]/F");
  jet_tree->Branch("chi2oh",_chi2oh,"chi2oh[24][64]/F");
  */
  if(_dotruthpar) jet_tree->Branch("truthparenergy",&_truthparenergy);
  if(_dotruthpar) jet_tree->Branch("truthpareta",&_truthpareta);
  if(_dotruthpar) jet_tree->Branch("truthparphi",&_truthparphi);
  if(_dotruthpar) jet_tree->Branch("truthparpt",&_truthparpt);
  if(_dotruthpar) jet_tree->Branch("truthparid",&_truthparid);
  
  if(_dowf)
    {
      _wft->Branch("runnum",&_runnum,"runnum/I");
      _wft->Branch("evtnum",&_evtnum,"evtnum/I");
      _wft->Branch("emwf",_emwf,"emwf[96][256][12]/i");
      _wft->Branch("ihwf",_ihwf,"ihwf[24][64][12]/i");
      _wft->Branch("ohwf",_ohwf,"ohwf[24][64][12]/i");
    }

  /*
  _wft->Branch("emieta",_emieta,"emieta[96][256]/I");
  _wft->Branch("ihieta",_ihieta,"ihieta[24][64]/I");
  _wft->Branch("ohieta",_ohieta,"ohieta[24][64]/I");
  _wft->Branch("emiphi",_emiphi,"emiphi[96][256]/I");
  _wft->Branch("ihiphi",_ihiphi,"ihiphi[24][64]/I");
  _wft->Branch("ohiphi",_ohiphi,"ohiphi[24][64]/I");

  jet_tree->Branch("emieta",_emieta,"emieta[96][256]/I");
  jet_tree->Branch("ihieta",_ihieta,"ihieta[24][64]/I");
  jet_tree->Branch("ohieta",_ohieta,"ohieta[24][64]/I");
  jet_tree->Branch("emiphi",_emiphi,"emiphi[96][256]/I");
  jet_tree->Branch("ihiphi",_ihiphi,"ihiphi[24][64]/I");
  jet_tree->Branch("ohiphi",_ohiphi,"ohiphi[24][64]/I");
  */
  if(_dowf)
    {
      _wft->Branch("emadcfit",_emadcfit,"emadcfit[96][256]/F");
      _wft->Branch("ihadcfit",_ihadcfit,"ihadcfit[24][64]/F");
      _wft->Branch("ohadcfit",_ohadcfit,"ohadcfit[24][64]/F");
      _wft->Branch("emt",_emt,"emt[96][256]/F");
      _wft->Branch("iht",_iht,"iht[24][64]/F");
      _wft->Branch("oht",_oht,"oht[24][64]/F");
    
      _wft->Branch("failscut",&_failscut,"failscut/I");
    }
  jet_tree->Branch("mbdavgt",_mbdavgt,"mbdavgt[2]/F");
  
  jet_tree->Branch("mbdhit",&_mbdhit,"mbdhit[2]/i");

  for(int i=0; i<nemx; ++i)
    {
      for(int j=0; j<nemy; ++j)
	{
	  _emieta[i][j] = i;
	  _emiphi[i][j] = j;
	  if(i< nx && j < ny)
	    {
	      _ihieta[i][j] = i;
	      _ihiphi[i][j] = j;
	      _ohieta[i][j] = i;
	      _ohiphi[i][j] = j;
	    }
	}
    }
  //jet_tree->Branch("nLayerEm",&_nLayerEm,"nLayerEm/I");
  //jet_tree->Branch("nLayerOh",&_nLayerOh,"nLayerOh/I");
  /*
  jet_tree->Branch("n2pc",&_n2pc,"n2pc/I");
  jet_tree->Branch("l2pcEta",&_l2pcEta,"l2pcEta/F");
  jet_tree->Branch("dPhi2pc",_dPhi2pc,"dPhi2pc[n2pc]/F");
  jet_tree->Branch("dEta2pc",_dEta2pc,"dEta2pc[n2pc]/F");
  jet_tree->Branch("dPhiLayer",_dPhiLayer,"dPhiLayer[jet_n]/F");
  */
  //jet_tree->Branch("emLayerJetPhi",_emLayerJetPhi,"emLayerJetPhi[nLayerEm]/F");
  //jet_tree->Branch("ohLayerJetPhi",_ohLayerJetPhi,"ohLayerJetPhi[nLayerOh]/F");
  //jet_tree->Branch("emLayerJetEta",_emLayerJetEta,"emLayerJetEta[nLayerEm]/F");
  //jet_tree->Branch("ohLayerJetEta",_ohLayerJetEta,"ohLayerJetEta[nLayerOh]/F");
  //jet_tree->Branch("emLayerJetET",_emLayerJetET,"emLayerJetET[nLayerEm]/F");
  //jet_tree->Branch("ohLayerJetET",_ohLayerJetET,"ohLayerJetET[nLayerOh]/F");
  
  //mbtree->Branch("mbevt",&_mbevt,"mbevt/I");
  
  _mbevt = 0;
  return Fun4AllReturnCodes::EVENT_OK;
}

//____________________________________________________________________________..
int zfinder::InitRun(PHCompositeNode *topNode)
{
  if(_debug > 1) cout << "Initializing!" << endl;
  return Fun4AllReturnCodes::EVENT_OK;
}

unsigned int get_towerindex(int packet, int pchan)
{
  return (192*(packet%1000-1))+pchan;
}

int get_bindex(int packet, int pchan, string what)
{
  unsigned int key = -1;
  if(packet/1000 == 6)
    {
      key = TowerInfoDefs::encode_emcal(get_towerindex(packet,pchan));
    }
  else if(packet/1000 == 7 || packet / 1000 == 8)
    {
      key = TowerInfoDefs::encode_hcal(get_towerindex(packet,pchan));
    }
  else
    {
      return -1;
    }

  if(what == "eta") return TowerInfoDefs::getCaloTowerEtaBin(key);
  else if(what == "phi") return TowerInfoDefs::getCaloTowerPhiBin(key);
  else return -1;
}

float zfinder::getEtaFromBin(int binEta)
{
  return ((2.2*binEta)/nx)-1.1;
}

float zfinder::getPhiFromBin(int binPhi)
{
  return ((2*M_PI*binPhi)/ny);}


float getEtaFromBinEM(int binEta)
{
  return ((2.2*binEta)/64)-1.1;
}

float getPhiFromBinEM(int binPhi)
{
  return ((2*M_PI*binPhi)/256);
}
//____________________________________________________________________________..

void print_debug(float jet_eta, float jet_phi, float tower_eta, float tower_phi, float dphi, float deta)
{
  cout << "printing debug info for dphi deta:" << endl;
  cout << "jet eta/phi: " << jet_eta << " " << jet_phi << endl;
  cout << "tower eta/phi: " << tower_eta << " " << tower_phi << endl;
  cout << "deta dphi:" << deta << " " << dphi << endl;
}

int zfinder::process_event(PHCompositeNode *topNode)
{

  

  if(_debug > 1) cout << endl << endl << endl << "zfinder: Beginning event processing" << endl;
  if(_nprocessed % 1000 == 0) cout << "processing event " << _nprocessed << endl;

  if(_isdat)
    {
      
      Gl1Packetv3* gl1 = findNode::getClass<Gl1Packetv3>(topNode, "14001");
      if(!gl1)
	{
	  cout << "No trigger info!" << endl;
	  return Fun4AllReturnCodes::ABORTRUN;
	}
      if(_debug > 1) cout << "zfinder: Getting gl1 trigger vector from: " << gl1 << endl;
      _triggervec = gl1->getScaledVector();
      
      int isjettrig = (_triggervec >> 22) & 1;
      int isjmbtrig = (_triggervec >> 18) & 1;
      //int ismbtrig = (_triggervec >> 10) & 1;
      //if(_debug > 2) cout << _triggervec << " " << isjettrig << " " << ismbtrig << endl;
      /*
	if(!isjettrig && !isjmbtrig)// && !ismbtrig)
	{
	if(_debug > 1) cout << "no jet trigger" << endl;
	return Fun4AllReturnCodes::ABORTEVENT;
	}
      */
      if(_nprocessed == 0)
	{
	  _prevraw18 = gl1->lValue(18,0);
	  _prevraw22 = gl1->lValue(22,0);
	  _prevlive18 = gl1->lValue(18,1);
	  _prevlive22 = gl1->lValue(22,1);
	  _isbadlive = 0;
	}
      long long unsigned int currentlive18, currentlive22, currentraw18, currentraw22;
      currentlive18 = gl1->lValue(18,1);
      currentlive22 = gl1->lValue(22,1);
      currentraw18 = gl1->lValue(18,0);
      currentraw22 = gl1->lValue(22,0);
      long long unsigned int live18diff = currentlive18 - _prevlive18;
      long long unsigned int live22diff = currentlive22 - _prevlive22;
      long long unsigned int raw18diff = currentraw18 - _prevraw18;
      long long unsigned int raw22diff = currentraw22 - _prevraw22;
      if(_nprocessed != 0)
	{
	  if(isjettrig && (((float)live22diff)/raw22diff < 0.1 || ((float)live18diff)/raw18diff < 0.1)) _isbadlive = 1;
	  else _isbadlive = 0;
	}
      else
	{
	  _isbadlive = 0;
	}
      
      _prevraw18 = currentraw18;
      _prevraw22 = currentraw22;
      _prevlive18 = currentlive18;
      _prevlive22 = currentlive22;
    }
  ++_nprocessed;
  /*
  if(ismbtrig)
    {
      if(_debug > 2) cout << "mb triggered" << endl;
      _mbevt++;
      if(!isjettrig) return Fun4AllReturnCodes::ABORTEVENT;
    }
  */
  if(_isdat)
    {
      int runnumber = 0;
      int evtnum = 0;
      EventHeader *runheader = findNode::getClass<EventHeader>(topNode, "EventHeader");
      if (!runheader)
	{
	  std::cout << "can't find runheader" << std::endl;
	  return Fun4AllReturnCodes::ABORTEVENT;
	}
      runnumber = runheader->get_RunNumber();
      evtnum = runheader->get_EvtSequence();
      _runnum = runnumber;
      _evtnum = evtnum;

    }
  PHNodeIterator itNode(topNode);
  PHCompositeNode* parNode = dynamic_cast<PHCompositeNode*>(itNode.findFirst("PHCompositeNode","PAR"));
  PdbParameterMap* flagNode;
  
  if(parNode) flagNode = findNode::getClass<PdbParameterMap>(parNode, "HasBeamBackground");
  else
    {
      cout << "No parNode! Abort run." << endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }
  
  
  if(flagNode) _cutParams.FillFrom(flagNode);
  else
    {
      cout << "No flagNode - abort run" << endl;
      return Fun4AllReturnCodes::ABORTRUN;
    }

  MbdVertexMap* mbdvtxmap = findNode::getClass<MbdVertexMapv1>(topNode, "MbdVertexMap");
  //GlobalVertexMap* gvtxmap = NULL; //findNode::getClass<GlobalVertexMapv1>(topNode, "GlobalVertexMap");

  float zvtx = NAN;
  /*
  if (gvtxmap)
    {
      if (gvtxmap->empty())
	{
	  if (_debug > 0)
	    {
	      std::cout << "gvtxmap empty - aborting event." << std::endl;
	    }
	  return Fun4AllReturnCodes::ABORTEVENT;
	}
      GlobalVertex *gvtx = gvtxmap->begin()->second;
      if (gvtx)
	{
	  auto startIter = gvtx->find_vertexes(_vtxtype);
	  auto endIter = gvtx->end_vertexes();
	  for (auto iter = startIter; iter != endIter; ++iter)
	    {
	      const auto &[type, vertexVec] = *iter;
	      if (type != _vtxtype)
		{
		  continue;
		}
	      for (const auto *vertex : vertexVec)
		{
		  if (!vertex)
		    {
		      continue;
		    }
		  zvtx = vertex->get_z();
		}
	    }
	}
      else
	{
	  if (_debug > 0)
	    {
	      std::cout << "gvtx is NULL! Aborting event." << std::endl;
	    }
	  return Fun4AllReturnCodes::ABORTEVENT;
	}
    }
  */
    if(mbdvtxmap)
    {
      for(auto iter = mbdvtxmap->begin(); iter != mbdvtxmap->end(); ++iter)
        {
          MbdVertex* mbdvtx = iter->second;
          if(mbdvtx) zvtx = mbdvtx->get_z();
          break;
        }
    }
    /*
  else
    {
      auto iter = gvtxmap->begin();
      while(iter != gvtxmap->end())
        {
          GlobalVertex* gvtx = iter->second;
          zvtx = gvtx->get_z();
          iter++;
          break;
        }
    }
  */
    /*
  if(std::isnan(zvtx) || abs(zvtx) > 30)
    {
      if(_debug > 1) cout << "no good zvtx!" << endl;
      return Fun4AllReturnCodes::ABORTEVENT;
    }
    */
  _zvtx = zvtx;
  if(std::isnan(_zvtx) || std::isnan(zvtx))
    {
      _zvtx = 0;
      zvtx = 0;
    }
  TowerInfoContainer *towers[3];
  towers[0] = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_CEMC_RETOWER");
  towers[1] = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_HCALIN");
  towers[2] = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_HCALOUT");
  TowerInfoContainer *regtows[3];
  regtows[0] = findNode::getClass<TowerInfoContainer>(topNode, "TOWERS_CEMC");
  regtows[1] = findNode::getClass<TowerInfoContainer>(topNode, "TOWERS_HCALIN");
  regtows[2] = findNode::getClass<TowerInfoContainer>(topNode, "TOWERS_HCALOUT");
  JetContainer *jets = findNode::getClass<JetContainerv1>(topNode, "AntiKt_Tower_HIRecoSeedsRaw_r04");//"AntiKt_unsubtracted_r04");
  if(!jets) jets = findNode::getClass<JetContainerv1>(topNode, "AntiKt_unsubtracted_r04");
  JetContainer* truthjets = NULL;
  if(!_isdat) truthjets = findNode::getClass<JetContainerv1>(topNode, "AntiKt_Truth_r04");
  //cout <<"truthjets: " <<  truthjets << endl;
  _tjet_n = 0;
  if(truthjets)
    {
      for(int i=0; i<truthjets->size(); ++i)
        {
          Jet* jet = truthjets->get_jet(i);
	  _tjet_pt[_tjet_n] = jet->get_pt();
	  if(_tjet_pt[_tjet_n] < 1) continue;
          _tjet_e[_tjet_n] = jet->get_e();
          _tjet_eta[_tjet_n] = jet->get_eta();
          _tjet_phi[_tjet_n] = jet->get_phi();
          _tjet_n++;
        }
    }

  
  if(_debug > 2) cout << towers[0] << " " << towers[1] << " " << towers[2] << endl;

  RawTowerGeomContainer *geom[3];
  geom[0] = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_CEMC");
  geom[1] = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALIN");
  geom[2] = findNode::getClass<RawTowerGeomContainer>(topNode, "TOWERGEOM_HCALOUT");


  MbdPmtContainer* mbdpmt = findNode::getClass<MbdPmtContainer>(topNode,"MbdPmtContainer");


  _mbdavgt[0] = 0;
  _mbdavgt[1] = 0;
  _mbdhit[0] = 0;
  _mbdhit[1] = 0;
  if(mbdpmt)
    {
      for(int i=0; i<128; ++i)
	{
	  MbdPmtHit* pmt = mbdpmt->get_pmt(i);
	  if(pmt)
	    {
	      if(pmt->get_q() > 0.4)
		{
		  ++_mbdhit[i/64];
		  _mbdavgt[i/64] += pmt->get_time();
		}
	    }
	}
      _mbdavgt[0]/=_mbdhit[0];
      _mbdavgt[1]/=_mbdhit[1];
    }
  else
    {
      if(_isdat)
	{
	  cout << "no MBD PMTs!!!" << endl;
	  return Fun4AllReturnCodes::ABORTRUN;
	}
    }

  







  
  //TowerInfoContainer* emcrt = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_CEMC_RETOWER");
  /*
  JetContainer* emjets = findNode::getClass<JetContainerv1>(topNode, "emtowjet");
  JetContainer* ohjets = findNode::getClass<JetContainerv1>(topNode, "ohtowjet");

  _nLayerEm = 0;
  _nLayerOh = 0;
  _n2pc = 0;

  if(emjets)
    {
      int tocheck = emjets->size();
      for(int i=0; i<tocheck; ++i)
	{
	  Jet* jet = emjets->get_jet(i);
	  if(jet)
	    {
	      float testJetET = jet->get_e()/cosh(jet->get_eta());
	      if(testJetET < 8) continue;
	      _emLayerJetEta[_nLayerEm] = jet->get_eta();
	      if(check_bad_jet_eta(_emLayerJetEta[_nLayerEm],zvtx,0.4)) continue;
	      _emLayerJetPhi[_nLayerEm] = jet->get_phi();
	      _emLayerJetET[_nLayerEm] = testJetET;
	      _nLayerEm++;
	    }
	}
    }

  if(ohjets)
    {
      int tocheck = ohjets->size();
      for(int i=0; i<tocheck; ++i)
	{
	  Jet* jet = ohjets->get_jet(i);
	  if(jet)
	    {
	      float testJetET = jet->get_e()/cosh(jet->get_eta());
	      if(testJetET < 5) continue;
	      _ohLayerJetEta[_nLayerOh] = jet->get_eta();
	      if(check_bad_jet_eta(_ohLayerJetEta[_nLayerOh],zvtx,0.4)) continue;
	      _ohLayerJetPhi[_nLayerOh] = jet->get_phi();
	      _ohLayerJetET[_nLayerOh] = testJetET;
	      _nLayerOh++;
	    }
	}
    }
  */
  /*
  int nchan = 1536;
  vector<vector<float>> emTowAbove1GeV;
  vector<vector<float>> ohTowAbove1GeV;
  _l2pcEta = 0;
  float maxTowET = 0;
  if(towers[0])
    {
      for(int i=0; i<nchan; ++i)
	{
	  TowerInfo* tower = towers[0]->get_tower_at_channel(i);
	  if(!tower->get_isGood()) continue;
	  int key = towers[0]->encode_key(i);
	  const RawTowerDefs::keytype geomkey = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::HCALIN, towers[0]->getTowerEtaBin(key), towers[0]->getTowerPhiBin(key));
	  RawTowerGeom *tower_geom = geom[1]->get_tower_geometry(geomkey);
	  float radius = 93.5;//tower_geom->get_center_radius();
	  float ihEta = tower_geom->get_eta();
	  float emZ = radius/(tan(2*atan(exp(-ihEta))));
	  float newz = emZ - zvtx;//tower_geom->get_center_z() - zvtx;
	  float newTheta = atan2(radius,newz);
	  float towerEta = -log(tan(0.5*newTheta));
	  float towerPhi = tower_geom->get_phi();
	  float towerET = tower->get_energy();
	  if(towerET < 1) continue;
	  if(towerET > maxTowET)
	    {
	      maxTowET = towerET;
	      _l2pcEta = towerEta;
	    }
	  vector<float> toPush = {towerEta, towerPhi};
	  emTowAbove1GeV.push_back(toPush);
	}
    }
  nchan = 1536;
  if(towers[2])
    {
      for(int i=0; i<nchan; ++i)
	{
	  TowerInfo* tower = towers[2]->get_tower_at_channel(i);
	  if(!tower->get_isGood()) continue;
	  int key = towers[2]->encode_key(i);
	  const RawTowerDefs::keytype geomkey = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::HCALOUT, towers[2]->getTowerEtaBin(key), towers[2]->getTowerPhiBin(key));
	  RawTowerGeom *tower_geom = geom[2]->get_tower_geometry(geomkey); //encode tower geometry                                                                                      

	  float radius = tower_geom->get_center_radius();
	  float newz = tower_geom->get_center_z() - zvtx;
	  float newTheta = atan2(radius,newz);
	  float towerEta = -log(tan(0.5*newTheta));
	  float towerPhi = tower_geom->get_phi();
	  //if(!_printedPhi) cout << "Phi tower : " << towerPhi << endl;
	  float towerET = tower->get_energy();
	  if(towerET < 1) continue;
	  if(towerET > maxTowET)
	    {
	      maxTowET = towerET;
	      _l2pcEta = towerEta;
	    }
	  vector<float> toPush = {towerEta, towerPhi};
	  ohTowAbove1GeV.push_back(toPush);
	}
    }
  _printedPhi = true;
  _n2pc = 0;
  for(int i=0; i<emTowAbove1GeV.size(); ++i)
    {
      for(int j=0; j<ohTowAbove1GeV.size(); ++j)
	{
	  if(_n2pc > 1000) break;
	  _dPhi2pc[_n2pc] = emTowAbove1GeV.at(i).at(1) - ohTowAbove1GeV.at(j).at(1);
	  if(_dPhi2pc[_n2pc] > M_PI) _dPhi2pc[_n2pc] -= 2*M_PI;
	  if(_dPhi2pc[_n2pc] < -M_PI) _dPhi2pc[_n2pc] += 2*M_PI;
	  _dEta2pc[_n2pc] = emTowAbove1GeV.at(i).at(0) - ohTowAbove1GeV.at(j).at(0);
	  ++_n2pc;
	}
    }
  */
  //float fracEM = 0;
  //float fracOH = 0;
  _maxTowE = 0;
  _subTowE = 0;
  _maxTowDiff = 0;
  _nBadChi2 = 0;
  float maxJetE = 0;
  float maxJetEta = 0;
  float maxJetPhi = 0;
  float subJetE = 0;
  float subJetPhi = 0;
  float maxLayerE[2] = {0};
  float maxEoverTot = 0;
  float SumEdEtadPhi = 0;
  float SumEdEtadPhiDiff = 0;
  float v1 = 0;
  float v2 = 0;
  float theta = 0;
  int jet_n = 0;
  _jet_n = 0;
  
  float jet_ecc = -1;
  float jet_lfrac = -1;
  for(int i=0; i<nx; ++i)
    {
      for(int j=0; j<ny; ++j)
	{
	  _jconem[i][j] = 0;
	  _jconih[i][j] = 0;
	  _jconoh[i][j] = 0;
	}
    }

  for(int i=0; i<100; ++i)
    {
      _jet_t[i] = 0;
      _jet_t_em[i] = 0;
      _jet_t_ih[i] = 0;
      _jet_t_oh[i] = 0;
    }
  //float Etot = 0;
  //bool isPerimeter[nx][ny] = {0};
  //float calE[3][nx][ny] = {0};
  //MinimumBiasInfov1* mb3 = findNode::getClass<MinimumBiasInfov1>(topNode,"mbc_bkgd3");
  //int maxCalN = -1;
  //float maxCaleE = 0;
  //float eTot = 0;

  float eccentricity = 0;
  //int failsceut = 0;
  if(jets)
    {
      int tocheck = jets->size();
      if(_debug > 2) cout << "Found " << tocheck << " jets to check..." << endl;
      for(int i=0; i<tocheck; ++i)
        {
          Jet *jet = jets->get_jet(i);
          if(jet)
            {
	      //float Etot = 0;
	      //if(_debug > 5) cout << "getting jet E/eta" << endl;
	      float testJetE = jet->get_pt();
	      float testJetPhi = jet->get_phi();
	      float sigEtaEta = 0;
	      float sigEtaPhi = 0;
	      float sigPhiPhi = 0;
	      //if(_debug > 5) cout << "jet E/eta: " << testJetE  << " " << jet->get_eta() << endl;
	      if(jet->get_pt() < 1) continue;
	      //if(_debug > 3) cout << "got a candidate jet" << endl;
	      _alljetfrcem[_jet_n] = 0;
	      _alljetfrcoh[_jet_n] = 0;
	      _jet_et[_jet_n] = jet->get_e();
	      _jet_eta[_jet_n] = jet->get_eta();
	      //if(check_bad_jet_eta(_jet_eta[_jet_n],zvtx,0.4)) continue;
	      _jet_pt[_jet_n] = jet->get_pt();
	      //if(abs(jet_eta[jet_n]) > 0.9) continue;
	      _jet_phi[_jet_n] = testJetPhi;//(jet->get_phi()>0?jet->get_phi():jet->get_phi()+2*M_PI);
	      TLorentzVector emAxis;
	      TLorentzVector ohAxis;
	      int ncomp = 0;
	      bool newMaxJetET = false;
	      if(testJetE > subJetE && testJetE < maxJetE)
		{
		  subJetE = testJetE;
		  subJetPhi = testJetPhi;
		}
	      if(_debug > 2) cout << "found a good jet!" << endl;
	      if(testJetE > maxJetE)
		{
		  newMaxJetET = true;
		  if(maxJetE > subJetE)
		    {
		      subJetE = maxJetE;
		      subJetPhi = maxJetPhi;
		    }
		  _maxTowE = 0;
		  _subTowE = 0;
		  _maxTowDiff = 0;
		  _nBadChi2 = 0;
		  _maxETowChi2 = 0;
		  _maxETowIsZS = -1;
		  _maxETowChi2Det = -1;
		  eccentricity = 0;
		  maxLayerE[0] = 0;
		  maxLayerE[1] = 0;
		  SumEdEtadPhi = 0;
		  SumEdEtadPhiDiff = 0;
		  v1 = 0;
		  v2 = 0;
		  maxJetE = testJetE;
		  maxJetEta = jet->get_eta();
		  maxJetPhi = jet->get_phi();//(jet->get_phi()>0?jet->get_phi():jet->get_phi()+2*M_PI);
		}
	
	      //if(_debug > 3) cout << "getting comp vec" << endl;
	      
	      for(int i=0; i<3; ++i)
		{
		  _maxTowChi2[i] = -1;
		  /*
		    for(int k=0; k<512; ++k)
		    {
		    _jetcompE[i][k] = 0;
		    _jetcompEta[i][k] = 0;
		    _jetcompPhi[i][k] = 0;
		    }
		  */
		}
	      float jet_t_Esum = 0;
	      float jet_emsum = 0;
	      float jet_ihsum = 0;
	      float jet_ohsum = 0;
	      for(auto comp: jet->get_comp_vec())
		{
		  ++ncomp;
		  unsigned int channel = comp.second;
		  TowerInfo* tower;
		  //cout << "type: " << comp.first << endl;
		  if(comp.first == 5 || comp.first == 26)
		    {
		      tower = towers[1]->get_tower_at_channel(channel);
		      float towerE = tower->get_energy();
		      float chi2 = tower->get_chi2();
		      //if(abs(chi2 - 0.08) < 0.00001 && towerE > 8) failscut = 2;
		      if(chi2 > _maxTowChi2[1]) _maxTowChi2[1] = chi2;
		      if(tower->get_isBadChi2()) _nBadChi2++;
		      //Etot += towerE;
		      int key = towers[1]->encode_key(channel);
		      const RawTowerDefs::keytype geomkey = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::HCALIN, towers[1]->getTowerEtaBin(key), towers[1]->getTowerPhiBin(key));
		      _jconih[towers[1]->getTowerEtaBin(key)][towers[1]->getTowerPhiBin(key)] = testJetE;
		      //if(_debug > 6) cout << "encoding tower geom" << endl;
		      RawTowerGeom *tower_geom = geom[1]->get_tower_geometry(geomkey); //encode tower geometry
		      float radius = tower_geom->get_center_radius();
		      float newz = tower_geom->get_center_z() - zvtx;
		      float newTheta = atan2(radius,newz);
		      float towerEta = -log(tan(0.5*newTheta));
		      float towerPhi = tower_geom->get_phi();

		      if(towerE > 0.1)
			{
			  jet_t_Esum += towerE;
			  jet_ihsum += towerE;
			  _jet_t_ih[_jet_n] += towerE*tower->get_time();
			  _jet_t[_jet_n] += towerE*tower->get_time();
			}
		      
		      /*
			_jetcompE[1][subcomp[1]] = towerE;
			_jetcompEta[1][subcomp[1]] = towerPhi;
			_jetcompPhi[1][subcomp[1]] = towerEta;
			subcomp[1]++;
		      */
		      //towerE /= cosh(towerEta);
		      if(towerE > _maxTowE)
			{
			  _subTowE = _maxTowE;
			  _maxTowE = towerE;
			  _maxETowChi2Det = 1;
			  _maxETowIsZS = tower->get_isZS();
			}
		      if(towerE < 0) continue;
		      //towerE = sqrt(towerE);
		      float dPhi = towerPhi - maxJetPhi;
		      if(dPhi > M_PI) dPhi -= 2*M_PI;
		      if(dPhi < -M_PI) dPhi += 2*M_PI;
		      float dEta = towerEta - maxJetEta;
		      //if(_debug > 4) cout << "IHCal" << endl;
		      //if(_debug > 4) print_debug(maxJetEta, maxJetPhi, towerEta, towerPhi, dPhi, dEta);
		      if(newMaxJetET)
			{
			  sigEtaEta += towerE*dEta*dEta;
			  sigEtaPhi += towerE*dEta*dPhi;
			  sigPhiPhi += towerE*dPhi*dPhi;
			  SumEdEtadPhi += 2*towerE*dPhi*dEta;
			  SumEdEtadPhiDiff += towerE*(dPhi*dPhi-dEta*dEta);
			}
		    }
		  else if(comp.first == 7 || comp.first == 27)
		    {
		      tower = towers[2]->get_tower_at_channel(channel);
		      float towerE = tower->get_energy();
		      
		      float chi2 = tower->get_chi2();
		      //if(abs(chi2 - 0.08) < 0.00001 && towerE > 8) failscut = 2;
		      if(chi2 > _maxTowChi2[2]) _maxTowChi2[2] = chi2;
		      if(tower->get_isBadChi2()) _nBadChi2++;
		      //fracOH += towerE;
		      //Etot += towerE;
		      int key = towers[2]->encode_key(channel);
		      const RawTowerDefs::keytype geomkey = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::HCALOUT, towers[2]->getTowerEtaBin(key), towers[2]->getTowerPhiBin(key));
		      _jconoh[towers[2]->getTowerEtaBin(key)][towers[2]->getTowerPhiBin(key)] = testJetE;
		      RawTowerGeom *tower_geom = geom[2]->get_tower_geometry(geomkey); //encode tower geometry
		      
		      float radius = tower_geom->get_center_radius();
		      float newz = tower_geom->get_center_z() - zvtx;
		      float newTheta = atan2(radius,newz);
		      float towerEta = -log(tan(0.5*newTheta));
		      float towerPhi = tower_geom->get_phi();

		      if(towerE > 0.1)
			{
			  jet_ohsum += towerE;
			  _jet_t_oh[_jet_n] += towerE*tower->get_time();
			  jet_t_Esum += towerE;
			  _jet_t[_jet_n] += towerE*tower->get_time();
			}
		      /*
			_jetcompE[2][subcomp[2]] = towerE;
			_jetcompEta[2][subcomp[2]] = towerPhi;
			_jetcompPhi[2][subcomp[2]] = towerEta;
			subcomp[2]++;
		      */
		      //towerE /= cosh(towerEta);
		      TLorentzVector tempOH;
		      tempOH.SetPtEtaPhiE(towerE,towerEta,towerPhi,tower->get_energy());
		      ohAxis += tempOH;
		      _alljetfrcoh[_jet_n] += towerE;
		      if(towerE > _maxTowE)
			{
			  _subTowE = _maxTowE;
			  _maxTowE = towerE;
			  _maxETowChi2Det = 2;
			  _maxETowIsZS = tower->get_isZS();
			}
		      if(newMaxJetET) maxLayerE[1] += towerE;
		      if(towerE < 0) continue;
		      //towerE = sqrt(towerE);
		      float dPhi = towerPhi - maxJetPhi;
		      if(dPhi > M_PI) dPhi -= 2*M_PI;
		      if(dPhi < -M_PI) dPhi += 2*M_PI;
		      float dEta = towerEta - maxJetEta;
		      //if(_debug > 4) cout << "OHCal" << endl;
		      //if(_debug > 4) print_debug(maxJetEta, maxJetPhi, towerEta, towerPhi, dPhi, dEta);
		      if(newMaxJetET)
			{
			  sigEtaEta += towerE*dEta*dEta;
			  sigEtaPhi += towerE*dEta*dPhi;
			  sigPhiPhi += towerE*dPhi*dPhi;
			  SumEdEtadPhi += 2*towerE*dPhi*dEta;
			  SumEdEtadPhiDiff += towerE*(dPhi*dPhi-dEta*dEta);
			}
		    }
		  else if(comp.first == 13 || comp.first == 28 || comp.first == 25)
		    {
		      tower = towers[0]->get_tower_at_channel(channel);
		      float towerE = tower->get_energy();
		      float chi2 = tower->get_chi2();
		      //if(abs(chi2 - 0.08) < 0.00001 && towerE > 8) failscut = 2;
		      if(chi2 > _maxTowChi2[0]) _maxTowChi2[0] = chi2;
		      if(tower->get_isBadChi2()) _nBadChi2++;
		      //Etot += towerE;
		      //fracEM += towerE;
		      int key = towers[0]->encode_key(channel);
		      const RawTowerDefs::keytype geomkey = RawTowerDefs::encode_towerid(RawTowerDefs::CalorimeterId::HCALIN, towers[0]->getTowerEtaBin(key), towers[0]->getTowerPhiBin(key));
		      _jconem[towers[0]->getTowerEtaBin(key)][towers[0]->getTowerPhiBin(key)] = testJetE;
		      RawTowerGeom *tower_geom = geom[1]->get_tower_geometry(geomkey); //encode tower geometry
		      
		      float radius = 93.5;//tower_geom->get_center_radius();
		      float ihEta = tower_geom->get_eta();
		      float emZ = radius/(tan(2*atan(exp(-ihEta))));
		      float newz = emZ - zvtx;//tower_geom->get_center_z() - zvtx;
		      float newTheta = atan2(radius,newz);
		      float towerEta = -log(tan(0.5*newTheta));
		      float towerPhi = tower_geom->get_phi();

		      if(towerE > 0.1)
			{
			  jet_emsum += towerE;
			  _jet_t_em[_jet_n] += towerE*tower->get_time();
			  jet_t_Esum += towerE;
			  _jet_t[_jet_n] += towerE*tower->get_time();
			}
		      
		      /*
			_jetcompE[0][subcomp[0]] = towerE;
			_jetcompEta[0][subcomp[0]] = towerPhi;
			_jetcompPhi[0][subcomp[0]] = towerEta;
			subcomp[0]++;
		      */
		      //towerE /= cosh(towerEta);
		      TLorentzVector tempEM;
		      tempEM.SetPtEtaPhiE(towerE,towerEta,towerPhi,tower->get_energy());
		      emAxis += tempEM;
		      _alljetfrcem[_jet_n] += towerE;
		      if(towerE > _maxTowE)
			{
			  _subTowE = _maxTowE;
			  _maxTowE = towerE;
			  _maxETowChi2Det = 0;
			  _maxETowIsZS = tower->get_isZS();
			}
		      if(newMaxJetET) maxLayerE[0] += towerE;
		      if(towerE < 0) continue;
		      //towerE = sqrt(towerE);
		      float dPhi = towerPhi - maxJetPhi;
		      if(dPhi > M_PI) dPhi -= 2*M_PI;
		      if(dPhi < -M_PI) dPhi += 2*M_PI;
		      float dEta = towerEta - maxJetEta;
		      //if(_debug > 4) cout << "EMCal" << endl;
		      //if(_debug > 4) print_debug(maxJetEta, maxJetPhi, towerEta, towerPhi, dPhi, dEta);
		      if(newMaxJetET)
			{
			  sigEtaEta += towerE*dEta*dEta;
			  sigEtaPhi += towerE*dEta*dPhi;
			  sigPhiPhi += towerE*dPhi*dPhi;
			  SumEdEtadPhi += 2*towerE*dPhi*dEta;
			  SumEdEtadPhiDiff += towerE*(dPhi*dPhi-dEta*dEta);
			}
		    }
		}
	      if(newMaxJetET)
		{
		  sigEtaEta /= maxJetE;
		  sigEtaPhi /= maxJetE;
		  sigPhiPhi /= maxJetE;
		  float a = sigEtaEta;
		  float b = sigEtaPhi;
		  float c = sigEtaPhi;
		  float d = sigPhiPhi;
		  float lam1 = (a+d+sqrt((a+d)*(a+d)-4*(a*d-b*c)))/2;
		  float lam2 = (a+d-sqrt((a+d)*(a+d)-4*(a*d-b*c)))/2;
		  v1 /= ncomp;
		  v2 /= ncomp;
		  maxEoverTot = (maxLayerE[0] > maxLayerE[1]?maxLayerE[0]:maxLayerE[1])/maxJetE;
		
		  //if(_debug > 3 && (lam1 < 0 || lam2 < 0)) cout << "lam1 or lam2 < 0, printing:" << lam1 << " " << lam2 << endl;
		  if(lam1 > lam2) eccentricity = 1-lam2/lam1;
		  else eccentricity = 1-lam1/lam2;
		  //if(_debug > 3) cout << "ecc/layer: " << eccentricity << " " << maxEoverTot << endl;
		}
	      //if(_debug > 3) cout << "ecc entries: " << h2_ecc_layer->GetEntries() << endl;
	      _alljetfrcoh[_jet_n] /= _jet_et[_jet_n];
	      _alljetfrcem[_jet_n] /= _jet_et[_jet_n];
	      _dPhiLayer[_jet_n] = emAxis.Phi() - ohAxis.Phi();
	      if(_dPhiLayer[_jet_n] > M_PI) _dPhiLayer[_jet_n] -= 2*M_PI;
	      if(_dPhiLayer[_jet_n] < -M_PI) _dPhiLayer[_jet_n] += 2*M_PI;
	      _jet_t[_jet_n] /= jet_t_Esum;
	      _jet_t_em[_jet_n] /= jet_emsum;
	      _jet_t_ih[_jet_n] /= jet_ihsum;
	      _jet_t_oh[_jet_n] /= jet_ohsum;
	      ++_jet_n;
	      if(_jet_n > 98) break;
	    }
	  else
	    {
	      continue;
	    }
	}
      if(!_jet_n)
	{
	  return Fun4AllReturnCodes::ABORTEVENT;
	}
      int failscut = 0;
      _bbfqavec = 0;
      _elmbgvec = 0;

      _bbfqavec = _cutParams.get_int_param("HasBeamBackground_StreakSidebandFilter") << 5;
      //if(_bbfqavec) cout << "bbfqavec: " << _bbfqavec <<  " and >> 5: " << (_bbfqavec >> 5) << endl;
	  //}
      _maxTowDiff = _maxTowE - _subTowE;
      //fracEM /= Etot;
      //fracOH /= Etot;
      float dphi = abs(maxJetPhi-subJetPhi);
      if(dphi > M_PI) dphi = 2*M_PI - dphi;
      //if(subJetE < 8) hfillnum=2;
      if(subJetE > 4) _isdijet = 1;
      else _isdijet = 0;
      /*
      if(failscut && _isdijet) _isdijet = 3;
      else if(failscut) _isdijet = 2;
      */
      //_isdijet += failscut;
      //{
      if(_debug > 3) cout << "assigning vars now" << endl;
      _eccentricity = eccentricity;
      _theta = theta;
      float maxEnergy = 0;
      float maxpt = 0;
      float subEnergy = 0;
      float subpt = 0;
      for(int k=0; k<_jet_n; ++k)
	{
	  if(_jet_pt[k] > maxpt)
	    {
	      subpt = maxpt;
	      subEnergy = maxEnergy;
	      maxpt = _jet_pt[k];
	      maxEnergy = _jet_et[k];
	    }
	  else if(_jet_pt[k] > subpt)
	    {
	      subpt = _jet_pt[k];
	      subEnergy = _jet_et[k];
	    }
	}
	      
      _frcoh = maxLayerE[1]/maxEnergy;
      _frcem = maxLayerE[0]/maxEnergy;
      _eta = maxJetEta;
      _phi = maxJetPhi;
      _jet_ET = maxEnergy;
      _dphi = dphi;
      _subjet_ET = subEnergy;
      if(_debug > 7) cout << "vars assigned" << endl;
      /*
      if(maxJetE > 4)
	{
	  if(_debug > 0) cout << "filling jet tree" << endl;
	  jet_tree->Fill();
	}
      */
      //}
      jet_ecc = eccentricity;
      jet_lfrac = maxEoverTot;
      failscut = (_bbfqavec >> 5) & 1;
      bool dPhiCut = (_dphi < 3*M_PI/4 || !_isdijet || maxJetE*0.2 > subJetE);
      bool loETCut = _frcem > 0.9 || _frcem < 0.1 || _frcoh < 0.1 || _frcoh > 0.9 || (1.-_frcem-_frcoh) > 0.5;// ((_frcem < 0.1) && (_jet_ET > (50*_frcem+20))) && (dPhiCut || !_isdijet);
      bool hiETCut = ((_frcem > 0.9) && (_jet_ET > (-50*_frcem+70))) && (dPhiCut || !_isdijet);
      bool ihCut = (_frcem+_frcoh) < 0.65;
      //bool fullCut = loETCut || hiETCut || ihCut || failscut || (_bbfqavec >> 5 & 0x1);
      int failsall = -1;//fullCut?1:0;
      if(!loETCut && !dPhiCut)
	{
	  failsall = 2;
	}
      else if(!loETCut)
	{
	  failsall = 1;
	}
      else if(!dPhiCut)
	{
	  failsall = 0;
	}
      if(_debug > 7) cout << "cuts set" << endl;



      

      if(_dotruthpar)
	{
	  _truthparenergy.clear();
	  _truthpareta.clear();
	  _truthparphi.clear();
	  _truthparpt.clear();
	  _truthparid.clear();
	  PHG4TruthInfoContainer *truthinfo = findNode::getClass<PHG4TruthInfoContainer>(topNode, "G4TruthInfo");
	  if(!truthinfo)
	    {
	      std::cout
		<< "MyJetAnalysis::process_event - Error can not find DST Truth Info node "
		<< "G4TruthInfo" << std::endl;
	      return Fun4AllReturnCodes::ABORTEVENT;
	    }
	  // Get the primary particle range
	  PHG4TruthInfoContainer::Range range = truthinfo->GetPrimaryParticleRange();

	  // Loop over the G4 truth (stable) particles
	  for (PHG4TruthInfoContainer::ConstIterator iter = range.first;
	       iter != range.second;
	       ++iter)
	    {
	      // Get this truth particle
	      const PHG4Particle *truth = iter->second;
	      
	      /// Get this particles momentum, etc.
	      float m_truthpx = truth->get_px();
	      float m_truthpy = truth->get_py();
	      float m_truthpz = truth->get_pz();
	      float m_truthpt_temp = sqrt(m_truthpx * m_truthpx + m_truthpy * m_truthpy);
	      if(m_truthpt_temp < 0.5) continue;
	      // float m_truthp = sqrt(m_truthpx * m_truthpx + m_truthpy * m_truthpy + m_truthpz * m_truthpz);
	      float m_tempe = truth->get_e();

	      _truthparenergy.push_back(m_tempe);
	      _truthparphi.push_back(atan(m_truthpy / m_truthpx));
	      _truthparpt.push_back(m_truthpt_temp);
	      float m_tempeta= atanh(m_truthpz / m_tempe);
	      /// Check for nans
	      if (!std::isfinite(m_tempeta))
		{
		  m_tempeta = -99;
		}
	      _truthpareta.push_back(m_tempeta);
	      _truthparid.push_back(truth->get_pid());	      
	    }
	}
      
      
      if((maxJetE > 10 && (!loETCut || !dPhiCut)) || maxJetE > 10 || _doall60 || !_isdat)
	{

	  for(int j=0; j<12; ++j)
	    {
	      for(int k=0; k<nemx; ++k)
		{
		  for(int l=0; l<nemy; ++l)
		    {
		      _emwf[k][l][j] = 0;
		    }
		}
	      for(int k=0; k<nx; ++k)
		{
		  for(int l=0; l<ny; ++l)
		    {
		      _ihwf[k][l][j] = 0;
		      _ohwf[k][l][j] = 0;
		    }
		}
	    }
	  
	  towers[0] = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_CEMC");
	  if(_debug > 1) cout << "adding to tree. towers: " << towers[0] << endl;
	  for(int j=0; j<3; ++j)
	    {
	      for(int k=0; k<(j==0?nemt:nt); ++k)
		{
		  TowerInfo* tower = towers[j]->get_tower_at_channel(k);
		  TowerInfo* regtow;
		  if(regtows[j]) regtow = regtows[j]->get_tower_at_channel(k);
		  else regtow=NULL;
		  //if(_debug > 9) cout << "tower: " << tower << endl;
		  int key = towers[j]->encode_key(k);
		  int eta = towers[j]->getTowerEtaBin(key);
		  int phi = towers[j]->getTowerPhiBin(key);
		  //if(_debug > 9) cout << "got key, eta, phi" << key << " " << eta << " " << phi << endl;
		  bool isbad = tower->get_isBadChi2();
		  bool isnocal = tower->get_isNoCalib();
		  bool ishot = tower->get_isHot();
		  if(j==0)
		    {
		      _emt[eta][phi] = tower->get_time();
		      _emtow[eta][phi] = tower->get_energy();
		      if(regtow) _emadcfit[eta][phi] = regtow->get_energy();
		      else _emadcfit[eta][phi] = -1;
		      _isbadem[eta][phi] = (isbad?1:0);
		      _ishotem[eta][phi] = (ishot?1:0);
		      _nocalem[eta][phi] = (isnocal?1:0);
		      _chi2em[eta][phi] = tower->get_chi2();
		    }
		  else if(j==1)
		    {
		      _iht[eta][phi] = tower->get_time();
		      _ihtow[eta][phi] = tower->get_energy();
		      if(regtow) _ihadcfit[eta][phi] = regtow->get_energy();
		      else _ihadcfit[eta][phi] = -1;
		      _isbadih[eta][phi] = (isbad?1:0);
		      _ishotih[eta][phi] = (ishot?1:0);
		      _nocalih[eta][phi] = (isnocal?1:0);
		      _chi2ih[eta][phi] = tower->get_chi2();
		    }
		  else if(j==2)
		    {
		      _oht[eta][phi] = tower->get_time();
		      _ohtow[eta][phi] = tower->get_energy();
		      if(regtow) _ohadcfit[eta][phi] = regtow->get_energy();
		      else _ohadcfit[eta][phi] = -1;
		      _isbadoh[eta][phi] = (isbad?1:0);
		      _ishotoh[eta][phi] = (ishot?1:0);
		      _nocaloh[eta][phi] = (isnocal?1:0);
		      _chi2oh[eta][phi] = tower->get_chi2();
		    }
		}
	      if(_debug > 3) cout << "did towers["<<j<<"]" << endl;
	    }
	  _failscut = failsall;
	  //drawCalo(towers, _jet_pt, _jet_eta, _jet_phi, _jet_n, jet_ecc, jet_lfrac, geom, zvtx, failsall, runnumber, evtnum, _frcoh, _frcem, maxLayerE[0], maxLayerE[1], maxJetE);
	  cout << "jet_tree = " << jet_tree << endl;
	  if(_dowf)
	    {
	      for(int p=6001; p<6129; ++p)
		{
		  if(_debug > 6) cout << "entered packet " << p << endl;
		  CaloPacketv1* packet = findNode::getClass<CaloPacketv1>(topNode,to_string(p));
		  if(!packet) continue;
		  for(int c=0; c<packet->getNrChannels(); ++c)
		    {		      
		      int etabin = get_bindex(p,c,"eta");
		      int phibin = get_bindex(p,c,"phi");
		      int ns = packet->getNrSamples();
		      if(_debug > 10)cout << "entered packet " << p << " channel " << c << " eta/phi " << etabin << " " << phibin << endl;
		      if(ns==12)
			{
			  for(int s=0; s<ns; ++s)
			    {
			      if(_debug > 12) cout << s << " ";
			      if(s>11)
				{
				  cout << "too many samples in " << p << " " << c  << endl;
				}
			      _emwf[etabin][phibin][s] = packet->getSample(c,s);
			    }
			}
		      else
			{
			  for(int s=0; s<12; ++s)
			    {
			      if(_debug > 12) cout << s << " ";
			      _emwf[etabin][phibin][s] = 0;
			    }
			}
		      cout << endl;
		    }
		}
	      
	      for(int p=7001; p<7009; ++p)
		{
		  if(_debug > 6) cout << "entered packet " << p << endl;
		  CaloPacketv1* packet = findNode::getClass<CaloPacketv1>(topNode,to_string(p));
		  if(!packet) continue;
		  for(int c=0; c<packet->getNrChannels(); ++c)
		    {
		      
		      int etabin = get_bindex(p,c,"eta");
		      int phibin = get_bindex(p,c,"phi");
		      int ns = packet->getNrSamples();
		      if(_debug > 10) cout << "entered packet " << p << " channel " << c  << " eta/phi " << etabin << " " << phibin << endl;
		      if(ns==12)
			{	      
			  for(int s=0; s<ns; ++s)
			    {
			      if(_debug > 12) cout << s << " ";
			      if(s>11)
				{
				  cout << "too many samples in " << p << " " << c  << endl;
				}
			      _ihwf[etabin][phibin][s] = packet->getSample(c,s);
			    }
			}
		      else
			{
			  for(int s=0; s<12; ++s)
			    {
			      if(_debug > 12) cout << s << " ";
			      _ihwf[etabin][phibin][s] = 0;
			    }
			}
		      cout << endl;
		    }
		}
	      
	      
	      for(int p=8001; p<8009; ++p)
		{
		  if(_debug > 6) cout << "entered packet " << p << endl;
		  CaloPacketv1* packet = findNode::getClass<CaloPacketv1>(topNode,to_string(p));
		  if(!packet) continue;
		  for(int c=0; c<packet->getNrChannels(); ++c)
		    {		      
		      int etabin = get_bindex(p,c,"eta");
		      int phibin = get_bindex(p,c,"phi");
		      if(_debug > 10) cout << "entered packet " << p << " channel " << c << " eta/phi " << etabin << " " << phibin << endl;
		      int ns = packet->getNrSamples();
		      if(ns==12 && !packet->getSuppressed(c))
			{			
			  for(int s=0; s<ns; ++s)
			    {
			      if(_debug > 12) cout << s << " ";
			      if(s>11)
				{
				  cout << "too many samples in " << p << " " << c  << endl;
				}
			      _ohwf[etabin][phibin][s] = packet->getSample(c,s);
			    }
			}
		      else
			{
			  for(int s=0; s<12; ++s)
			    {
			      if(_debug > 12) cout << s << " ";
			      _ohwf[etabin][phibin][s] = 0;
			    }
			}
		      cout << endl;
		    }
		}
	      cout << "done writing to arrays" << endl;
	      _wft->Fill();
	    }
	  jet_tree->Fill();
	  cout << "drew calo" << endl;
	}
      
      
    }
  else
    {
      if(_debug > 0) cout << "no jets" << endl;
      return Fun4AllReturnCodes::ABORTEVENT;
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
  if (Verbosity() > 0)
    {
      std::cout << "zfinder::End(PHCompositeNode *topNode) This is the End..." << std::endl;
    }
  cout << "ending run" << endl;
  if(jet_tree->GetEntries() > 0)
    {
      if(!_doall60 && _isdat) _f = new TFile(_filename.c_str(), "RECREATE");
      cout << "file created" << endl;
      if(!_doall60 && _isdat) _f->cd();
      cout << "cded to file" << endl;
      if(!_doall60 && _isdat) jet_tree->SetDirectory(_f);
      cout << "tree set to directory of file" << endl;
      jet_tree->Write();
      cout << "tree written" << endl;
      _f->Write();
      cout << "file written" << endl;
      _f->Close();
      cout << "file closed" << endl;
    }
  if(_dowf)
    {
      if(_wft->GetEntries() > 0)
	{
	  if(!_doall60) _wff = new TFile(_wfilename.c_str(),"RECREATE");
	  if(!_doall60) _wff->cd();
	  if(!_doall60) _wft->SetDirectory(_wff);
	  _wft->Write();
	  _wff->Write();
	  _wff->Close();
	}
    }
  cout << "file saved if necessary" << endl;
  //delete jet_tree;
  //delete _f;
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
