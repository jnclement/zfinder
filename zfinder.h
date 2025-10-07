#ifndef ZFINDER_H
#define ZFINDER_H
#include <fun4all/SubsysReco.h>
#include <calobase/RawTowerContainer.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/TowerInfoContainer.h>
class PHCompositeNode;
class CentralityInfo;
class zfinder : public SubsysReco
{
 public:

  zfinder(const std::string &name = "zfinder", const int debug = 0);

  virtual ~zfinder();

  float new_eta(int channel, TowerInfoContainer* towers, RawTowerGeomContainer* geom, RawTowerDefs::CalorimeterId caloID, float testz);
  
  int Init(PHCompositeNode *topNode) override;

  int InitRun(PHCompositeNode *topNode) override;

  int process_event(PHCompositeNode *topNode) override;

  int ResetEvent(PHCompositeNode *topNode) override;

  int EndRun(const int runnumber) override;

  int End(PHCompositeNode *topNode) override;

  int Reset(PHCompositeNode * /*topNode*/) override;

  void Print(const std::string &what = "ALL") const override;

 
 private:
  int _debug;
  std::string _name;
  float _zvtx;
  unsigned int _nprocessed;
};

#endif // ZFINDER
