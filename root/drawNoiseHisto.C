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
  return pchs->get(icha) == 1;
}

struct NoiseSpecifier {
  string sspec;     // e.g. taiXX-50
  string sjob;      // e.g. taiXX
  string sreco;     // e.g. cor
  string snsam;     // e.g. 50
  NoiseSpecifier(string sspecIn) {
    sspec = sspecIn;
    string::size_type ipos1 = sspec.find("-");
    string srem = sspec;
    if ( ipos1 != string::npos ) {
      snsam = srem.substr(ipos1+1);
      srem = srem.substr(0, ipos1);
    }
    sjob = srem;
    sreco = srem.substr(0, 3);
  }
  bool isCoherent() const { return snsam.size(); }
  bool isKeScale() const {
    return sjob.find("adc") == string::npos;
  }
};

int fillFromTps(string jobsuf, const NoiseSpecifier& nspec, string sfrun, string sfevts, int itps, HistMap& hsts) {
  string myname = "fillFromTps: ";
  int iapas[6] = {3, 5, 2, 6, 1, 4};
  int iapa = iapas[itps];
  string sapa = std::to_string(iapa);
  string stps = std::to_string(itps);
  // Find the noise summary histograms.
  Index nvie = 0;
  NameVector  sfils;
  IndexVector nvchs;
  string sfilpre_noapa = nspec.sjob + jobsuf +                      "/run" + sfrun + "evts"
                         + sfevts + "/chmet_samRms" + nspec.snsam + "_";
  string sfilpre_byapa = nspec.sjob + jobsuf + "/dpcr_apa" + sapa + "/run" + sfrun + "evts"
                         + sfevts + "/chmet_samRms" + nspec.snsam + "_";
  string sfilsuf = "_run" + sfrun + ".tpad";
  string sviews[4] = {"u", "v", "z", "c"};
  // Make list of names to try.
  NameVector sfilcans;
  sfilcans.push_back(sfilpre_noapa + "tps" + stps + sfilsuf);
  sfilcans.push_back(sfilpre_byapa + "tps" + stps + sfilsuf);
  sfilcans.push_back(sfilpre_noapa + "apa" + sapa + sfilsuf);
  sfilcans.push_back(sfilpre_byapa + "apa" + sapa + sfilsuf);
  sfilcans.push_back(sfilpre_noapa + "tpp" + stps + sviews[0] + sfilsuf);
  sfilcans.push_back(sfilpre_byapa + "tpp" + stps + sviews[0] + sfilsuf);
  Index ifil;
  for ( ifil=0; ifil<sfilcans.size(); ++ifil ) {
    Name sfil = sfilcans[ifil];
    if ( gSystem->AccessPathName(sfil.c_str()) == 0 ) break;
  }
  if ( ifil < 4 ) {
    nvie = 1;
    nvchs.push_back(2560);
    sfils.push_back(sfilcans[ifil]);
  } else if ( ifil < 6 ) {
    nvie = 4;
    nvchs.push_back(800);
    nvchs.push_back(800);
    nvchs.push_back(480);
    nvchs.push_back(480);
    string sfilpre = ifil==4 ? sfilpre_noapa : sfilpre_byapa;
    for ( Index ivie=0; ivie<nvie; ++ivie ) sfils.push_back(sfilpre + "tpp" + stps + sviews[ivie] + sfilsuf);
  } else {
    cout << myname << "ERROR: Unable to find any of the following files:" << endl;
    for ( string badfil : sfilcans ) {
      cout << myname << "ERROR:   " << badfil << endl;
    }
    return 1;
  }
  cout << myname << "===== Processing " << nvie << " view" << (nvie == 1 ? "" : "s") << endl;
  for ( Index ivie=0; ivie<nvie; ++ivie ) {
    string sfil = sfils[ivie];
    cout << myname << "Processing " << sfil << endl;
    Index nvch = nvchs[ivie];
    const TPadManipulator* pmanin = TPadManipulator::read(sfil);
    if ( pmanin == nullptr ) {
      cout << myname << "ERROR: Unable to read pad file:" << endl;
      cout << myname << "  " << sfil << endl;
      return 1;
    }
    // Fetch graph.
    TGraph* pg = pmanin->graph();
    if ( pg->GetN() != nvch ) {
      cout << myname << "Graph has unexpected bin count: " << pg->GetN()
           << " !' " << nvch << endl;
      return 4;
    }
    for ( HistMap::value_type ihst : hsts ) {
      string hnam = ihst.first;
      TH1* phf = ihst.second;
      bool doBad = hnam.find("Bad") != string::npos;
      bool doGood = hnam.find("Good") != string::npos;
      Index ivch1 = 0;
      Index ivch2 = nvch;
      if ( hnam.substr(0,3) == "all" ) {
      } else if ( hnam.substr(0,2) == "uv" ) {
        if ( nvie == 4 ) {
          if ( ivie > 1 ) continue;
        } else if ( nvie == 1 ) {
          ivch2 = 1600;
        }
      } else if ( hnam.substr(0,2) == "zc" ) {
        if ( nvie == 4 ) {
          if ( ivie < 2 ) continue;
        } else if ( nvie == 1 ) {
          ivch1 = 1600;
        }
      } else {
        cout << myname << "Unexpected histo name: " << hnam << endl;
      }
      cout << myname << "Filling " << hnam << "(" << doGood << ", " << doBad << ")" << endl;
      bool noisy = 0;
      for ( Index ivch=ivch1; ivch<ivch2; ++ivch ) {
        double xch, val;
        if ( pg->GetPoint(ivch, xch, val) != int(ivch) ) {
          cout << "ERROR: Point not found: " << ivch << endl;
          continue;
        }
        Index ich = xch + 0.1;
        bool isbad = isBad(ich);
        if ( doGood && isbad ) continue;
        if ( doBad && !isbad ) continue;
        if ( noisy ) cout << ich << ": " << val << endl;
        phf->Fill(val);
      }
    }
  }
  return 0;
}

// Create a pad for one pair of histograms.
int plotNoiseHisto(string jobsuf, string sspec, float ymaxin, string sfrun, string sfevts, string dopt, TPadManipulator* pman) {
  string myname = "drawNoiseHisto: ";
  bool titleAsLabel = true;
  bool biglabs = true;
  NoiseSpecifier nspec(sspec);
  int itpss[6] = {0, 1, 2, 3, 4, 5};
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
  double xmax = 14.0;
  string sxttl = "Noise [(ADC count)/tick]";
  string sntype = "sample";
  if ( nspec.isKeScale() ) {
    if ( nspec.isCoherent() ) {
      xmax = 5.0;
      sxttl = "Noise [ke/(" + nspec.snsam + " tick)]";
      sntype = "coherent";
      sntype = "integrated";
      //nbin = 80;
    } else {
      xmax = 0.35;
      sxttl = "Noise [ke/tick]";
    }
  } else if ( nspec.isCoherent() ) {
    xmax = 200.0;
    sxttl = "Noise [(ADC count)/(" + nspec.snsam + " tick)]";
    sntype = "coherent";
    sntype = "integrated";
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
    int rstat = fillFromTps(jobsuf, nspec, sfrun, sfevts, itps, hsts);
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
  double dyl = 0.05;
  double yl2 = 0.85;
  double yl1 = yl2 - (hnams.size()+0.5)*dyl;
  TLegend* pleg = new TLegend(0.55, yl1, 1.00, yl2);
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
    ssmean << std::fixed << fac*ph->GetMean();
    string slab = legdescs[hnam] + ": #LT#sigma#GT = " + ssmean.str() + " " + sunit;
    pleg->AddEntry(ph, slab.c_str(), "l");
  }
  // Create labels.
  vector<TLatex*> labs;
  double xlab = 0.55;
  double ylab = 0.35;
  double dylab = 0.06;
  labs.push_back(new TLatex(xlab, ylab, "#bf{ProtoDUNE-SP}"));
  ylab -= dylab;
  if ( titleAsLabel ) {
    labs.push_back(new TLatex(xlab, ylab, sttlSum.c_str()));
    sttlSum = "";
    ylab -= dylab;
  }
  string slab = sfrun;
  while ( slab.size() && slab[0] == '0' ) slab = slab.substr(1);
  slab = "Run " + slab;
  labs.push_back(new TLatex(xlab, ylab, slab.c_str()));
  // Set margins and label sizes.
  if ( sttlSum.size() == 0 ) pman->setMarginTop(0.05);
  if ( biglabs ) {
    pman->setLabelSizeX(0.05);
    pman->setLabelSizeY(0.05);
    pman->setMarginLeft(0.13);
    pman->setMarginBottom(0.12);
  }

  // Draw.
  if ( doMan ) {
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
    if ( ymaxin > 0.0 ) ymax = ymaxin;
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
TPadManipulator* drawNoiseHisto(string jobsuf, string outpre, string a_sspec, string sfrun, string sfevts, string dopt, string printSufs =".png") {
  string myname = "drawNoiseHisto: ";
  StringManipulator sman(a_sspec, true);
  StringManipulator::StringVector sspecs = sman.split(",");
  int nssp = sspecs.size();
  if ( nssp == 0 ) return nullptr;
  TPadManipulator* pmantop = new TPadManipulator(1400, 1000);
  if ( nssp > 1 ) {
    int dim = int(sqrt(nssp-1)) + 1;
    pmantop->split(dim, dim);
  }
  string fnam = outpre + "noise_";
  for ( int issp=0; issp<nssp; ++issp ) {
    string sspec = sspecs[issp];
    cout << myname << sspec << endl;
    StringManipulator smflds(sspec, false);
    StringManipulator::StringVector sflds = smflds.split(":");
    sspec = sflds[0];
    float ymax = 0.0;
    if ( sflds.size() > 1 ) {
      istringstream ssymax(sflds[1]);
      ssymax >> ymax;
      cout << myname << "Set ymax to " << ymax << " for " << sspec << endl;
    }
    if ( issp ) fnam += "-";
    fnam += sspec;
    TPadManipulator* pman = pmantop->man(issp);
    int pstat = plotNoiseHisto(jobsuf, sspec, ymax, sfrun, sfevts, dopt, pman);
    if ( pstat ) {
      cout << myname << "Pad fill returned error " <<  pstat << endl;
      delete pmantop;
      return nullptr;
    }
  }
  if ( printSufs.size() ) {
    if ( dopt.size() ) fnam += "_" + dopt;
    fnam += "_run" + sfrun + "evts" + sfevts + ".{" + printSufs + "}";;
    cout << myname << "Printing " << fnam << endl;
    pmantop->print(fnam);
  }
  return pmantop;
}
