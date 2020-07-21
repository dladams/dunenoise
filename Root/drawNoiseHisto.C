// drawNoiseHisto.C
//
// David Adams
// June 2020
//
// Script to create summary noise distribution histograms.
//
// srec is expected to be format RRRQQ, RRRQQ-WW, RRRQQ:PP, RRRQQ-WW:PP
//   RRR = cor, cni, ...
//   QQ = qualifier, e.g. w1
//   WW = # samples 
//   PP = plane specifier = one or more of u, v, z, c
//
// The job directory is
//   RRRjobsub/run.../
// or
//   RRRjobsub/dpcr_CR/run.../

using HistMap = std::map<string, TH1*>;
using IntMap = std::map<string, int>;
using Name = std::string;
using NameVector = std::vector<Name>;
using NameMap = std::map<Name, Name>;
using Index = unsigned int;
using IndexVector = vector<Index>;

bool isBad(unsigned int icha) {
  const string myname = "isBad: ";
  DuneToolManager* ptm = DuneToolManager::instance();
  const IndexMapTool* pchs = nullptr;
  if ( ptm == nullptr ) {
    cout << myname << "ERROR: Tool manager not found." << endl;
    return true;
  } else {
    pchs = ptm->getShared<IndexMapTool>("channelStatus");
    if ( pchs == nullptr ) {
      cout << myname << "ERROR: Tool not found: channelStatus." << endl;
      return true;
    }
  }
  // 13-may-2020
  //return pchs->get(icha) == 1;
  return pchs->get(icha) > 0;
}

// sjob-svar:ymax
struct NoiseSpecifier {
  string sspec;     // e.g. tai-nsgms50, cnr-nsgrms:200
  string sreco;     // e.g. tai, cnr
  string svar;      // e.g. nsgrms, nsgrms50
  string snsam;     // e.g. "", "50"
  int nsam = 0;
  float ymax = 0.0;
  bool doPlot = true;
  NoiseSpecifier(string sspecIn) {
    const string myname = "NoiseSpecifier::ctor: ";
    sspec = sspecIn;
    string::size_type ipos1 = sspec.find("-");
    if ( ipos1 != string::npos ) {
      svar = sspec.substr(ipos1+1);
      sreco = sspec.substr(0, ipos1);
    }
    StringManipulator smflds(svar, false);
    StringManipulator::StringVector sflds = smflds.split(":");
    if ( sflds.size() > 1 ) {
      svar = sflds[0];
      if ( sflds[1] == "noplot" ) {
        doPlot = false;
      } else {
        istringstream ssymax(sflds[1]);
        ssymax >> ymax;
        cout << myname << "Set ymax to " << ymax << " for " << sspec << endl;
      }
    }
    if ( svar.substr(0, 6) == "nsgrms" && svar.size() > 6 ) {
      string snsam = svar.substr(6);
      istringstream ssnsam(snsam);
      ssnsam >> nsam;
    }
  }
  bool isKeScale() const {
    if ( svar.find("nsgrms") == string::npos ) return false;
    return sreco.find("adc") == string::npos;
  }
  bool isAdcScale() const {
    if ( svar.find("nsgrms") == string::npos ) return false;
    return sreco.find("adc") != string::npos;
  }
  bool isUtc() const {
    return svar.find("utc") != string::npos;
  }
  string recvar() const { return sreco + "-" + svar; }
};

int fillFromTps(string filpat, const NoiseSpecifier& nspec, string sfrun, string sdet, int itps, HistMap& hsts) {
  string myname = "fillFromTps: ";
  int iapas[6] = {3, 5, 2, 6, 1, 4};
  int iapa = iapas[itps];
  string sapa = std::to_string(iapa);
  string stps = std::to_string(itps);
  // Find the noise summary histograms.
  string::size_type ipos = sfrun.find("-");
  string sfrun1 = sfrun.substr(0, ipos);
  string sfrun2 = sfrun.substr(ipos+1);
  StringManipulator sman(filpat, true);
  sman.replace("%RUN%", sfrun);
  sman.replace("%RUN1%", sfrun1);
  sman.replace("%RUN2%", sfrun2);
  sman.replace("%REC%", nspec.sreco);
  sman.replace("%VAR%", nspec.svar);
  string infile = sman.str();
  string sviews[4] = {"u", "v", "z", "c"};
  if ( gSystem->AccessPathName(infile.c_str()) != 0 ) {
    cout << myname << "File not found: " << infile << endl;
    return 1;
  }
  int ichmod = 1;
  int ichcol = 0;
  if ( sdet == "pdsp" ) {
    ichmod = 2560;
    ichcol = 1600;
  } else if ( sdet == "iceberg3" ) {
    ichmod = 1280;
    ichcol =  800;
  } else {
    cout << myname << "Invalid detector: " << sdet << endl;
    return 3;
  }
  cout << myname << "===== Processing " << infile << endl;
  const TPadManipulator* pmanin = TPadManipulator::read(infile);
  if ( pmanin == nullptr ) {
    cout << myname << "ERROR: Unable to read pad file:" << endl;
    cout << myname << "  " << infile << endl;
    return 1;
  }
  // Fetch graph.
  TGraph* pg = pmanin->graph();
  int nvch = pg->GetN();
  for ( HistMap::value_type ihst : hsts ) {
    string hnam = ihst.first;
    TH1* phf = ihst.second;
    bool doBad = hnam.find("Bad") != string::npos;
    bool doGood = hnam.find("Good") != string::npos;
    Index ivch1 = 0;
    Index ivch2 = nvch;
    bool isInd = hnam.substr(0,2) == "uv";
    bool isCol = hnam.substr(0,2) == "zc";
    cout << myname << "Filling " << hnam << "(" << doGood << ", " << doBad << ")" << endl;
    bool noisy = 0;
    for ( Index ivch=ivch1; ivch<ivch2; ++ivch ) {
      double xch, val;
      if ( pg->GetPoint(ivch, xch, val) != int(ivch) ) {
        cout << "ERROR: Point not found: " << ivch << endl;
        continue;
      }
      Index ich = xch + 0.1;
      Index ichred = ich%ichmod;
      if ( isInd && ichred >= ichcol ) continue;
      if ( isCol && ichred <  ichcol ) continue;
      bool isbad = isBad(ich);
      if ( doGood && isbad ) continue;
      if ( doBad && !isbad ) continue;
      if ( noisy ) cout << ich << ": " << val << endl;
      phf->Fill(val);
    }
  }
  return 0;
}

// Create a pad for one pair of histograms.
int plotNoiseHisto(string filpat, string sspec, string sfrun, string sdet, string dopt,
                   TPadManipulator* pman, ostream& txtout) {
  string myname = "plotNoiseHisto: ";
  cout << myname << filpat << ", " << sspec << endl;
  bool titleAsLabel = true;
  bool biglabs = true;
  NoiseSpecifier nspec(sspec);
  vector<int> itpss;
  string explab;
  if ( sdet == "pdsp" ) {
    vector<int> tmp = {0, 1, 2, 3, 4, 5};
    itpss = tmp;
    explab = "#bf{ProtoDUNE-SP}";
  } else if ( sdet == "iceberg3" ) {
    itpss.push_back(0);
    explab = "#bf{Iceberg} run 3";
  } else {
    cout << myname << "Invalid detector: " << sdet << endl;
    return 3;
  }
  HistMap hsts;
  NameMap descs;
  NameMap legdescs;
  descs["all"] = "all";
  descs["uv"] = "uv";
  descs["zc"] = "zc";
  descs["uvGood"] = "good uv";
  descs["zcGood"] = "good zc";
  descs["allBad"] = "bad";
  legdescs["uvGood"] = "Induction";
  legdescs["zcGood"] = "Collection";
  IntMap cols;
  LineColors lc;
  cols["all"] = lc.black();
  cols["uv"] = lc.red();
  cols["uvGood"] = lc.red();
  cols["zc"] = lc.blue();
  cols["zcGood"] = lc.blue();
  cols["allBad"] = lc.gray();
  NameVector hnamsAll = {"all", "uv", "zc", "uvGood", "zcGood", "allBad"};
  NameVector hnams;
  string srem = dopt;
  while ( srem.size() ) {
    string::size_type ipos = srem.find("-");
    string hnam = srem.substr(0, ipos);
    bool isHnam = false;
    for ( string hnamChk : hnamsAll ) if ( hnam == hnamChk ) isHnam = true;
    if ( isHnam ) {
      hnams.push_back(hnam);
    } else {
      cout << myname << "Invalid channels specifier option: " << hnam << endl;
      return 1;
    }
    if ( ipos == string::npos ) break;
    srem = srem.substr(ipos+1);
    cout << "XXX: " << hnam << " , " << srem << " has length " << srem.size() << endl;
  }
  int nbin = 70;
  double xmin = 0.0;
  double xmax = 0.0;
  string sxttl = "Noise [(ADC count)/tick]";
  string sntype = "sample";
  if ( nspec.isKeScale() ) {
    if ( nspec.nsam ) {
      xmax = 5.0;
      sxttl = "Noise [ke/(" + nspec.snsam + " tick)]";
      sntype = "integrated";
    } else {
      xmax = 0.35;
      sxttl = "Noise [ke/tick]";
    }
  } else if ( nspec.isAdcScale() ) {
    if ( nspec.nsam ) {
      xmax = 14.0;
      sxttl = "Noise [(ADC count)/tick]";
      sntype = "integrated";
    } else {
      xmax = 200.0;
      sxttl = "Noise [(ADC count)/(" + nspec.snsam + " tick)]";
    }
  }
  string sttlSuf;
  if ( nspec.sreco == "cor" ) sttlSuf = " without CNR";
  if ( nspec.sreco == "tai" ) sttlSuf = " w/o CNR";
  if ( nspec.sreco == "cni" ) sttlSuf = " with CNR";
  if ( nspec.sreco == "cnr" ) sttlSuf = " with CNR";
  string sttlSum = sntype + " noise" + sttlSuf;
  sttlSum[0] = toupper(sttlSum[0]);
  int lineWidth = 2;
  for ( string hnam : hnams ) {
    string sttl = "TPC " + sntype + " noise for " + descs[hnam] + " channels" + sttlSuf;
    sttl += "; " + sxttl;
    sttl += "; # channels";
    TH1* ph = new TH1F(hnam.c_str(), sttl.c_str(), nbin, xmin, xmax);
    ph->SetDirectory(nullptr);
    ph->SetLineWidth(lineWidth);
    ph->SetLineColor(cols[hnam]);
    hsts[hnam] = ph;
  }
  for ( int itps : itpss ) {
    int rstat = fillFromTps(filpat, nspec, sfrun, sdet, itps, hsts);
    if ( rstat ) {
      cout << myname << "Exiting because fillFromTps returned " << rstat << endl;
      return 2;
    }
  }
  cout << myname << "Histogram count: " << hsts.size() << endl;
  for ( HistMap::value_type ent : hsts ) {
    TH1* ph = ent.second;
    cout << ph->GetName() << ": " << ph->GetEntries() << endl;
  }
  // Display and print plots.
  bool doMan = pman != nullptr;
  string hopt = "hist";
  // Add overflows to the last bin.
  bool showOver = true;
  if ( showOver ) {
    for ( string hnam : hnams ) {
      TH1* ph = hsts[hnam];
      int nbin = ph->GetNbinsX();
      double nover = ph->GetBinContent(nbin+1);
      if ( nover > 0.0 ) {
        double nsum = ph->GetBinContent(nbin) + nover;
        ph->SetBinContent(nbin, nsum);
      }
    }
  }
  // Create legend.
  bool showMedian = true;
  double dyl = 0.05;
  double yl2 = 0.85;
  double yl1 = yl2 - (hnams.size()+0.5)*dyl;
  double xl1 = showMedian ? 0.40 : 0.55;
  TLegend* pleg = new TLegend(xl1, yl1, 0.95, yl2);
  pleg->SetBorderSize(0);
  pleg->SetFillStyle(0);
  pleg->SetMargin(0.10);  // Fraction used for symbol
  int isty = 0;
  int prec = nspec.isKeScale() ? 0 : 1;
  double fac = nspec.isKeScale() ? 1000.0 : 1.0;
  string sunit = nspec.isKeScale() ? "e" : "ADC";
  for ( string hnam : hnams ) {
    TH1* ph = hsts[hnam];
    ph->SetLineStyle(++isty);
    ostringstream ssmean;
    ssmean.precision(prec);
    double mean = fac*ph->GetMean();
    ssmean << std::fixed << mean;
    string slab = legdescs[hnam] + ": #bar{#sigma} = " + ssmean.str() + " " + sunit;
    double median = 0.0;
    if ( showMedian ) {
      double q = 0.5;
      ph->GetQuantiles(1, &median, &q);
      median *= fac;
      ostringstream ssmed;
      ssmed.precision(prec);
      ssmed << std::fixed << median;
      slab += ", #tilde{#sigma} = " + ssmed.str() + " " + sunit;
    }
    pleg->AddEntry(ph, slab.c_str(), "l");
    txtout << setw(15) << nspec.recvar() << " " << setw(10) << hnam << " ";
    int wval = 12;
    if ( nspec.isUtc() ) {
      txtout << setw(wval) << int(mean+0.499);
      if ( showMedian ) txtout << " " << setw(wval) << int(median+0.499);
    } else {
      txtout << setw(wval) << mean;
      if ( showMedian ) txtout << " " << setw(wval) << median;
    }
    txtout << "\n";
  }
  // Create labels.
  vector<TLatex*> labs;
  double xlab = 0.55;
  double ylab = 0.45;
  double dylab = 0.06;
  labs.push_back(new TLatex(xlab, ylab, explab.c_str()));
  ylab -= dylab;
  if ( titleAsLabel ) {
    labs.push_back(new TLatex(xlab, ylab, sttlSum.c_str()));
    sttlSum = "";
    ylab -= dylab;
  }
  string slab = sfrun;
  while ( slab.size() && slab[0] == '0' ) slab = slab.substr(1);
  string::size_type ipos = slab.find("-");
  if ( ipos == string::npos ) {
    slab = "Run " + slab;
  } else {
    string srun1 = slab.substr(0, ipos);
    string srun2 = slab.substr(ipos+1);
    while ( srun2.size() && srun2[0] == '0' ) srun2 = srun2.substr(1);
    slab = "Runs " + srun1 + " - " + srun2;
  }
  labs.push_back(new TLatex(xlab, ylab, slab.c_str()));
  // Draw.
  if ( doMan ) {
    // Set margins and label sizes.
    if ( sttlSum.size() == 0 ) pman->setMarginTop(0.05);
    if ( biglabs ) {
      pman->setLabelSizeX(0.05);
      pman->setLabelSizeY(0.05);
      pman->setMarginLeft(0.13);
      pman->setMarginBottom(0.12);
    }
    double ymax = 0.0;
    for ( string hnam : hnams ) {
      TH1* ph = hsts[hnam];
      double ymaxh = ph->GetMaximum();
      if ( ymaxh > ymax ) ymax = ymaxh;
      ph->SetLineWidth(lineWidth);
      pman->add(ph, hopt);
      hopt = "hist same";
      lineWidth = 3;
    }
    ymax *= 1.03;
    if ( nspec.ymax > 0.0 ) ymax = nspec.ymax;
    for ( TLatex* ptxt : labs ) {
      ptxt->SetNDC();
      ptxt->SetTextFont(42);
      pman->add(ptxt);
    }
    pman->setTitle(sttlSum);
    pman->setRangeY(0.0, ymax);
    pman->add(pleg);
    pman->addAxis();
    pman->centerAxisTitles();
  } else {
    TCanvas* pcan = new TCanvas;
    for ( string hnam : hnams ) {
      TH1* ph = hsts[hnam];
      ph->Draw(hopt.c_str());
      hopt = "hist same";
    }
    pleg->Draw();
  }
  return 0;
}

// Create a pad with multiple histograms.
// filpat is the pattern for the input file
//        %RUN% is the run number or range from sfrun
//        %RUN1%, %RUN2% are the run numbers from sfrun=RUN1-RUN2
//        %REC% is the reco level from specs
//        %VAR% is the variable level from specs
// outnam is the base for out put files with the following substitutions
//        %RUN% is the run number or range from sfrun
//        %ALLSPEC% is built from the specifier from each plot
// a_sspec - comma-separated list of specifiers: REC-VAR, REC-N:YMAX
//           REC = tai, cnr
//           VAR = nsgrms, nsgrms50, ...
//           YMAX  = limit for plot
//           One plot is filled for each entry.
// sfrun is string rep of run number as it appears in file name
// sdet is pdsp or iceberg3
// scurves is which curves to draw for each plot, e.g. "zcGood-uvGood"
// plotSufs is the tpad suffix for output plots
TPadManipulator* drawNoiseHisto(string filpat, string outnam, string a_sspec, string sfrun,
                                string sdet, string scurves, string plotSufs =".png") {
  string myname = "drawNoiseHisto: ";
  StringManipulator sman(a_sspec, true);
  StringManipulator::StringVector sspecs = sman.split(",");
  int nssp = sspecs.size();
  int nman = 0;
  if ( nssp == 0 ) return nullptr;
  // Build allspec.
  string allspec;
  for ( int issp=0; issp<nssp; ++issp ) {
    string sspec = sspecs[issp];
    NoiseSpecifier nspec(sspec);
    if ( issp ) allspec += "_";
    allspec += nspec.recvar();  // sspec without ymax
    if ( nspec.doPlot ) ++nman;
  }
  cout << "nssp, nman = " << nssp << ", " << nman << endl;
  // Build the output file base.
  StringManipulator smano(outnam, true);
  smano.replace("%RUN%", sfrun);
  smano.replace("%CURVES%", scurves);
  smano.replace("%ALLSPEC%", allspec);
  string fnam = smano.str();
  // Build the manipulator.
  TPadManipulator* pmantop = new TPadManipulator(1400, 1000);
  if ( nman > 1 ) {
    int dim = int(sqrt(nman-1) + 0.01) + 1;
    pmantop->split(dim, dim);
  }
  string txtnam = fnam + ".txt";
  ofstream txtout(txtnam.c_str());
  int iman = 0;
  for ( int issp=0; issp<nssp; ++issp ) {
    string sspec = sspecs[issp];
    cout << myname << sspec << endl;
    NoiseSpecifier nspec(sspec);
    TPadManipulator* pman = nspec.doPlot ? pmantop->man(iman++) : nullptr;
    int pstat = plotNoiseHisto(filpat, sspec, sfrun, sdet, scurves, pman, txtout);
    if ( pstat ) {
      cout << myname << "Pad fill returned error " <<  pstat << endl;
      delete pmantop;
      return nullptr;
    }
  }
  cout << "Output summary: " << txtnam << endl;
  if ( plotSufs.size() ) {
    fnam += ".{" + plotSufs + "}";
    cout << myname << "Printing " << fnam << endl;
    pmantop->print(fnam);
  }
  return pmantop;
}
