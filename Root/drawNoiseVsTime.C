// drawNoiseVsTime.C
//
// David Adams
// May 2020
//
// Draw noise vs time. and after CNR plots.

TPadManipulator* drawNoiseVsTime(int irun1, int irun2, string spatin, string sxvar, string snoise, float ymax) {
  const string myname = "drawNoiseVsTime: ";
  TGraph* pgr = new TGraph();
  map<string, TGraph*> grs;
  LineColors lc;
  string sruns;
  TGraph* pgrmax = new TGraph;
  pgrmax->SetMarkerColor(0);
  pgrmax->SetMarkerStyle(0);
  pgrmax->GetXaxis()->SetTitle("Run number");
  pgrmax->GetYaxis()->SetTitle("Noise [e]");
  TLegend* pleg = new TLegend(0.80, 0.15, 0.95, 0.35);
  pleg->SetBorderSize(0);
  pleg->SetFillStyle(0);
  pleg->SetMargin(0.10);  // Fraction used for symbol
  for ( int irun=irun1; irun<=irun2; ++irun ) {
    StringManipulator sman(spatin, true);
    string srun = to_string(irun);
    while ( srun.size() < 6 ) srun = "0" + srun;
    if ( irun == irun1 ) sruns = srun;
    if ( irun == irun2 ) sruns += "-" + srun;
    sman.replace("%RUN%", srun);
    string sfin = sman.str();
    ifstream fin(sfin);
    if ( ! fin ) {
      cout << myname << "Unable to open " << sfin << endl;
      continue;
    }
    cout << myname << "Reading " << sfin << endl;
    string svar;
    string scrv;
    double val;
    float valmax = 0.0;
    while ( fin >> svar >> scrv >> val ) {
      string::size_type ipos = svar.find("-");
      string svar1 = svar.substr(0, ipos);
      string svar2 = svar.substr(ipos + 1);
      if ( svar2 == snoise ) {
        string gname = svar + ":" + scrv;
        if ( grs.count(gname) == 0 ) {
          grs[gname] = new TGraph;
          TGraph* pgr = grs[gname];
          int icol = lc.black();
          string slab = scrv;
          if ( slab.substr(0,2) == "uv" ) slab = "Induction";
          if ( slab.substr(0,2) == "zc" ) slab = "Collection";
          if ( svar1 == "tai" ) slab += " w/o CNR";
          if ( svar1 == "cnr" ) slab += " with CNR";
          if ( scrv[0] == 'z' || scrv['0'] == 'c' ) icol = lc.blue();
          if ( scrv[0] == 'u' || scrv['0'] == 'v' ) icol = lc.red();
          int imrk = 2;
          if ( svar1 == "cnr" ) imrk = 27;
          pgr->SetMarkerColor(icol);
          pgr->SetMarkerStyle(imrk);
          pleg->AddEntry(pgr, slab.c_str(), "p");
        }
        TGraph* pgr = grs[gname];
        pgr->SetPoint(pgr->GetN(), irun, val);
        if ( val > valmax ) valmax = val;
        cout << myname << svar1 << " " << svar2 << " " << scrv << ": " << val << endl;
      }
      if ( valmax > 0.0 ) pgrmax->SetPoint(pgrmax->GetN(), irun, valmax);
    }
    //vector<string> svar = "tai-nsgrms50";
  }
  string sttl;
  if ( snoise == "nsgrms" ) sttl = "Sample noise";
  else if ( snoise.substr(0,6) == "nsgrms" ) sttl = "Integrated noise";
  TPadManipulator* ppad = new TPadManipulator(1400, 500);
  ppad->add(pgrmax, "P");
  for ( auto& igr : grs ) {
    TGraph* pgr = igr.second;
    ppad->add(pgr, "P");
    ppad->addAxis();
    ppad->setGridY();
  }
  ppad->add(pleg);
  string fnamout = "noiseVsTime_" + snoise + "_" + sruns + ".png";
  ppad->setRangeX(irun1-1, irun2+1);
  ppad->setRangeY(0, ymax);
  ppad->showGraphOverflow();
  ppad->setTitle(sttl);
  //ppad->centerAxisTitles();
  ppad->print(fnamout);
  return ppad;
}
