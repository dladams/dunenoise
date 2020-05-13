//  ApaWirePlot_tool.cc
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

#include "ApaWirePlot.h"
#include "dune/Geometry/IcebergChannelGeo.h"
#include <iostream>

using std::cout;
using std::endl;
using std::string;

//**********************************************************************

ApaWirePlot::ApaWirePlot(const fhicl::ParameterSet& ps)
: m_LogLevel(ps.get<Index>("LogLevel")),
  m_Xwidth(ps.get<Index>("Xwidth")),
  m_Ywidth(ps.get<Index>("Ywidth")),
  m_Xmin(ps.get<float>("Xmin")),
  m_Xmax(ps.get<float>("Xmax")),
  m_Zmin(ps.get<float>("Zmin")),
  m_Zmax(ps.get<float>("Zmax")),
  m_Ymin(ps.get<float>("Ymin")),
  m_Ymax(ps.get<float>("Ymax")),
  m_pad(m_Xwidth, m_Ywidth) {
  const string myname = "ApaWirePlot::ctor: ";
  TH1* ph = new TH2F("hax", ";z [cm];y [cm]", 1, m_Zmin, m_Zmax, 1, m_Ymin, m_Ymax);
  ph->SetStats(0);
  m_pad.add(ph, "h");
  cout << myname << "  LogLevel: " << m_LogLevel << endl;
  cout << myname << "      Xmin: " << m_Xmin << endl;
  cout << myname << "      Xmax: " << m_Xmax << endl;
  cout << myname << "      Zmin: " << m_Zmin << endl;
  cout << myname << "      Zmax: " << m_Zmax << endl;
  cout << myname << "      Ymin: " << m_Ymin << endl;
  cout << myname << "      Ymax: " << m_Ymax << endl;
}

//**********************************************************************

int ApaWirePlot::addChannel(Index icha, int icol, int isty, int iwid) {
  IcebergChannelGeo cg(icha);
  int iwire = 0;
  for ( IcebergChannelGeo::EndPoints eps : cg.wires() ) {
    float x1 = eps.first.x();
    string sstat = "Skip";
    float y1 = eps.first.y();
    float z1 = eps.first.z();
    float x2 = eps.first.x();
    float y2 = eps.second.y();
    float z2 = eps.second.z();
    if ( x1 > m_Xmin && x1 < m_Xmax ) {
      TLine lin(z1, y1, z2,y2);
      lin.SetLineColor(icol);
      lin.SetLineWidth(iwid);
      m_pad.add(&lin);
      sstat = "Keep";
    }
    if ( m_LogLevel >= 3 ||
         (m_LogLevel == 2 && sstat == "Keep") ) {
      cout << sstat << " " << icha << "-" << iwire
           <<     ": (" << x1 << ", " << y1 << ", " << z1
           << ") --> (" << x2 << ", " << y2 << ", " << z2
           << endl;
    }
    ++iwire;
  }
  return 0;
}
  

//**********************************************************************

int ApaWirePlot::print(string fnam, string ttl) {
  m_pad.setTitle(ttl);
  m_pad.addAxis();
  return m_pad.print(fnam);
}

//**********************************************************************
