#ifndef CALOANA_H__
#define CALOANA_H__

#include <fun4all/SubsysReco.h>

// Forward declarations
class Fun4AllHistoManager;
class PHCompositeNode;
class TFile;
class TNtuple;

class CaloAna : public SubsysReco
{
 public:
  //! constructor
  CaloAna(const std::string &name = "CaloAna", const std::string &fname = "MyNtuple.root");

  //! destructor
  virtual ~CaloAna();

  //! full initialization
  int Init(PHCompositeNode *);

  //! event processing method
  int process_event(PHCompositeNode *);

  //! end of run method
  int End(PHCompositeNode *);

  int process_g4hits(PHCompositeNode *);
  int process_g4cells(PHCompositeNode *);
  int process_towers(PHCompositeNode *);
  int process_clusters(PHCompositeNode *);

  void Detector1(const std::string &name) { detector1 = name; }
  void Detector2(const std::string &name) { detector2 = name; }
  void Detector3(const std::string &name) { detector3 = name; }
  void Detector4(const std::string &name) { detector4 = name; }

 protected:
  std::string detector1;
  std::string detector2;
  std::string detector3;
  std::string detector4;
  
  std::string outfilename;
  Fun4AllHistoManager *hm;
  TFile *outfile;
  TNtuple *g4hitntuple;
  TNtuple *g4cellntuple;
  TNtuple *towerntuple;
  TNtuple *clusterntuple;
  int evtid;
};

#endif
