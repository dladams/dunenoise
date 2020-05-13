// ibPlotWires.C
//
// David Adams
// May 2020
//
// Plot the wires for one plane on one side of the Iceberg APA.
//
// view: Wire orientations: u, v or z
// side: front is z< 0, back is z > 0
// fnam: Name for output file e.g. det.{png,tpad}
// ttl = Plot title
// wids = width for each channel; empty uses 4 for all

void ibPlotWires(string view, string side, string fnam, string ttl,
                 vector<int> wids = vector<int>(0)) {
  LineColors lc;
  vector<int> fcols(11,0);
  for ( int icol=0; icol<10; ++icol ) fcols[icol+1] = lc.color(icol, 10);
  string tname = "ibWirePlot";
  if ( side == "front" ) tname += "Front";
  else if ( side == "back" ) tname += "Back";
  else if ( side == "all" ) tname += "All";
  else {
    cout << "Invalid side: " << side << endl;
    return;
  }
  std::unique_ptr<ApaWirePlot> pwp = ptm->getPrivate<ApaWirePlot>(tname);
  int chanoff = 0;
  int chsign = 1;
  vector<int> ifmbs;
  int nfch = 40;
  if ( view == "u" ) {
    chanoff = 0;
    chsign = -1;
    for ( int ifmb= 5; ifmb>=1; --ifmb ) ifmbs.push_back(ifmb);
    for ( int ifmb=10; ifmb>=6; --ifmb ) ifmbs.push_back(ifmb);
  } else if ( view == "v" ) {
    chanoff = 400;
    for ( int ifmb=1; ifmb<=10; ++ifmb ) ifmbs.push_back(ifmb);
  } else if ( view == "z" ) {
    nfch = 48;
    if ( side == "front" ) {
      chanoff = 800;
      for ( int ifmb=1; ifmb<=5; ++ifmb ) ifmbs.push_back(ifmb);
    } else if ( side == "back" ) {
      chanoff = 1040;
      chsign = -1;
      for ( int ifmb=10; ifmb>=6; --ifmb ) ifmbs.push_back(ifmb);
    }
  } else {
    cout << "Invalid view: " << view << endl;
    return;
  }
  int icha = chanoff;
  for ( int ifmb : ifmbs ) {
    int icol = fcols[ifmb];
    for ( int ifch=0; ifch<nfch; ++ifch ) {
      int wid = 4;
      if ( wids.size() ) {
        wid = icha < wids.size() ? wids[icha] : 0.0;
      }
      if ( wid > 0.0 ) pwp->addChannel(icha, icol, 1, wid);
      ++icha;
    }
  }
  int wsave = gStyle->GetLineWidth();
  int ssave = gStyle->GetLineScalePS();
  gStyle->SetLineWidth(5);
  gStyle->SetLineScalePS(1);
  pwp->print(fnam, ttl);
  gStyle->SetLineWidth(wsave);
  gStyle->SetLineScalePS(ssave);
}
