#ifndef CALOVTXRECO_H
#define CALOVTXRECO_H
#include <fun4all/SubsysReco.h>
#include <calobase/RawTowerGeomContainer.h>
#include <calobase/TowerInfoContainer.h>
#include <globalvertex/CaloVertexMapv1.h>
class PHCompositeNode;
class CaloVtxReco : public SubsysReco
{
 public:

  CaloVtxReco(const std::string &name = "CaloVtxReco", const int debug = 0);

  virtual ~CaloVtxReco();

  int createNodes(PHCompositeNode *topNode);

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
  CaloVertexMap* _calovtxmap;
};

#endif // CALOVTXRECO
