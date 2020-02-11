bool widescale() { return  false; }

// jobdir - Input files are in PROC{jobdir}/dftpowlog...
// outpre - Ouput plots are {outpre}dftmerge...
// sprcin - Comma-separated List of PROC values, e.g. "tai,cnr"
// srun - string representation of run number, e.g. "5240"
// slabin - labe prefixes for the PROCs, e.g. "After tail removal:,After CNR:"
// ylogmin - Min value for y-axis in log plot, e.g. 1.e-5
// ylogmax - Max value for y-axis in log plot, e.g. 0.004
// ylinmax - Max value for y-axis in linear plot, e.g. 0.0024
void mergeDftPlanes(string jobdir, string outpre,
                    string sprcin, string srun, string slabin,
                    double ylogmin, double ylogmax, double ylinmax) {
  string myname = "mergeplanes: ";
  bool useescale = ylinmax > 1.0;
  bool biglabs = true;
  bool showtitle = false;
  StringManipulator smprcs(sprcin, false);
  StringManipulator::StringVector sprcs = smprcs.split(",");
  int nprc = sprcs.size();
  cout << myname << "Process count is " << nprc << endl;
  //if ( slabin.size() == 0 ) slabin = "After tail removal:,After CNR:";
  if ( slabin.size() == 0 ) slabin = "W/o CNR:,With CNR:";
  StringManipulator smlabs(slabin, false);
  StringManipulator::StringVector slabs = smlabs.split(",");
  int nlab = slabs.size();
  cout << myname << "Process label count is " << nlab << endl;
  TPadManipulator omantop(1400, 1000);
  omantop.split(2,2);
  double xlmin = 0.58;
  double xlmax = 0.93;
  double ylmax = 0.88;
  double ylmin = ylmax - 0.07*(nprc+0.5);
  if ( ylmin < 0.40 ) ylmin = 0.40;
  string sprclab;
  //double ymax[4] = {0.002, 0.008, 0.15, 0.008};
  double ymax[4] = {0.0024, 0.0024, 0.0024, 0.0024};
  if ( ylinmax > 0.0 ) for ( double& ymaxval : ymax ) ymaxval = ylinmax;
  vector<string> padlabs = {"C wires", "U wires", "Z wires", "V wires"};
  for ( int iprc=0; iprc<nprc; ++iprc ) {
    string sprc = sprcs[iprc];
    StringManipulator smflds(sprc, false);
    string dopt = "HIST";
    int icol = LineColors::color(iprc);
    if ( iprc ) dopt += " SAME";
    if ( iprc ) sprclab += "-";
    sprclab += sprc;
    string fnam = sprc + jobdir + "/dftpowtlog_run00" + srun + ".tpad";
    TPadManipulator* pmantop = TPadManipulator::read(fnam);
    if ( pmantop == nullptr ) {
      cout << myname << "Unable to find tpm in " << fnam << endl;
      return;
    }
    for ( int iman=0; iman<4; ++iman ) {
      TPadManipulator* pman = pmantop->man(iman);
      TH1* ph = (TH1*) pman->hist()->Clone();
      if ( useescale ) {
        ph->GetYaxis()->SetTitle("Power/tick/channel [e^{2}/(20 kHz)]");
      } else {
        ph->GetYaxis()->SetTitle("Power/tick/channel [(ke)^{2}/(20 kHz)]");
      }
      int nobj = pman->objects().size();
      int nhst = 1;
      string sttl = pman->getTitle();
      for ( int iobj=0; iobj<nobj; ++iobj ) {
        TH1* pho = pman->getHist(iobj);
        if ( pho == nullptr ) continue;
        ph->Add(pho);
        ++nhst;
      }
      float sfac = 1.0/nhst;
      if ( useescale ) sfac *= 1.e6;
      ph->Scale(sfac);
      ph->SetLineColor(icol);
      ph->Print();
      if ( nprc == 2 && iprc == 0 ) {
        //ph->SetLineStyle(2);
        ph->SetLineWidth(3);
      }
      TPadManipulator* poman = omantop.man(iman);
      poman->centerAxisTitles();
      poman->add(ph, dopt);
      TH1* pho = dynamic_cast<TH1*>(poman->lastObject());
      if ( iprc == 0 ) {
        poman->addLegend(xlmin, ylmin, xlmax, ylmax);
      }
      cout << myname << "Extra object count for pad " << iman << ": " << poman->objects().size() << endl;
      float ytxt = 0.42;
      if ( widescale() ) ytxt = 0.36;
      //float ytxt = 0.36;
      if ( iprc + 1 == nprc ) {
        poman->setLogY();
        poman->setRangeY(0.0, ymax[iman]);
        poman->setLogRangeY(ylogmin, ylogmax);
        poman->addAxis();
        //poman->setGrid();
        if ( showtitle ) {
          poman->setTitle("DFT power");
        } else {
          poman->setTitle("");
          poman->setMarginTop(0.05);
        }
        double xlab = 0.65;
        if ( widescale() ) xlab = 0.70;
        TLatex* ptxt = new TLatex(xlab, ytxt, "#bf{ProtoDUNE-SP}");
        ptxt->SetNDC();
        ptxt->SetTextFont(42);
        poman->add(ptxt);
        ytxt -= 0.06;
        ptxt = new TLatex(xlab, ytxt, padlabs[iman].c_str());
        ptxt->SetNDC();
        ptxt->SetTextFont(42);
        poman->add(ptxt);
        ytxt -= 0.06;
        string stxt = "Run " + srun;
        ptxt = new TLatex(xlab, ytxt, stxt.c_str());
        ptxt->SetNDC();
        ptxt->SetTextFont(42);
        poman->add(ptxt);
      }
      // Set margins.
      if ( biglabs ) {
        poman->setLabelSizeX(0.05);
        poman->setLabelSizeY(0.05);
        poman->setMarginLeft(0.14);
        poman->setMarginBottom(0.12);
      }
      TLegend* pleg = poman->getLegend();
      pleg->SetMargin(0.12);  // Fraction used for symbol
      if ( pleg == nullptr ) {
        cout << myname << "ERROR: Legend not found." << endl;
      } else {
        string slab = iprc < nlab ? slabs[iprc] : sprc;
        double powsum = 0.0;
        for ( int ibin=0; ibin<=pho->GetNbinsX(); ++ibin ) {
          powsum += pho->GetBinContent(ibin);
        }
        double powrms = sqrt(powsum);
        ostringstream sspow;
        sspow.precision(0);
        if ( ! useescale ) powrms *= 1000;
        sspow << "#sqrt{#Sigma} = " << fixed << powrms << " e";
        slab += " " + sspow.str();
        pleg->AddEntry(pho, slab.c_str(), "l");
      }
    }
  }
  string ssuf = ".{png,pdf,tpad}";
  string sofnam = outpre + "dftmergelog_" + sprclab + "_run" + srun + ssuf;
  cout << myname << "Printing " << sofnam << endl;
  omantop.print(sofnam);
  for ( int iman=0; iman<4; ++iman ) {
    TPadManipulator* pman = omantop.man(iman);
    pman->setLogY(false);
    if ( ! biglabs ) pman->setMarginLeft(0.12);
    //TH1* phframe = pman->frameHist();
    //cout << "Frame X title: " << phframe->GetXaxis()->GetTitle() << endl;
  }
  sofnam = outpre + "dftmergelin_" + sprclab + "_run" + srun + ssuf;
  omantop.update();
  cout << myname << "Printing " << sofnam << endl;
  omantop.print(sofnam);
}
