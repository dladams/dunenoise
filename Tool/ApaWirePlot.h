//  ApaWirePlot.h
//
// David Adams
// May 2019
//
// Draw wires for an APA in the plane parallel to the APA.
//
// Parameters:
//   LogLevel: Noise level: 0=nonee, ...
//   Xmin, Xmax: Wires in the range are drawn
//   Zmin, Zmax: Horizontal limits for the plot
//   Ymin, Ymax: Vertical limits for the plot

#include "fhiclcpp/ParameterSet.h"
#include "art/Utilities/ToolMacros.h"
#include "dune/DuneCommon/TPadManipulator.h"
#include "TH2F.h"
#include <string>

class ApaWirePlotInterface {
public:
  using Index = unsigned int;
  virtual ~ApaWirePlotInterface() =default;
  virtual int addChannel(Index icha, int icol=1, int isty=1, int iwid=1) =0;
  virtual int print(std::string fname, std::string ttl) =0;
  virtual TPadManipulator& pad() =0;
};

class ApaWirePlot : public ApaWirePlotInterface {
public:
  explicit ApaWirePlot(fhicl::ParameterSet const& ps);
  int addChannel(Index icha, int icol=1, int isty=1, int iwid=1) override;
  int print(std::string fnam, std::string ttl) override;
  TPadManipulator& pad() override { return m_pad; };
private:
  // Configuration parameters.
  Index m_LogLevel;
  Index m_Xwidth;
  Index m_Ywidth;
  float m_Xmin;
  float m_Xmax;
  float m_Zmin;
  float m_Zmax;
  float m_Ymin;
  float m_Ymax;
  // Derived data.
  TPadManipulator m_pad;
};

DEFINE_ART_CLASS_TOOL(ApaWirePlot)

