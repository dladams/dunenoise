bool widescale() { return  false; }

// jobdir - Input files are in PROC{jobdir}/dftpowlog...
// outpre - Ouput plots are {outpre}dftmerge...
// sprcin - Comma-separated List of PROC values, e.g. "tai,cnr"
// srun - string representation of run number, e.g. "5240"
// slabin - labe prefixes for the PROCs, e.g. "After tail removal:,After CNR:"
// ylogmin - Min value for y-axis in log plot, e.g. 1.e-5
// ylogmax - Max value for y-axis in log plot, e.g. 0.004
void mergeDftPlanes(string jobdir, string outpre,
                    string sprcin, string srun, string slabin,
                    double ylogmin, double ylogmax) {
  string myname = "mergeplanes: ";
  StringManipulator smprcs(sprcin, false);
  StringManipulator::StringVector sprcs = smprcs.split(",");
  int nprc = sprcs.size();
  cout << myname << "Process count is " << nprc << endl;
  if ( slabin.size() == 0 ) slabin = "After tail removal:,After CNR:";
  StringManipulator smlabs(slabin, false);
  StringManipulator::StringVector slabs = smlabs.split(",");
  int nlab = slabs.size();
  cout << myname << "Process label count is " << nlab << endl;
  TPadManipulator omantop(1400, 1000);
  omantop.split(2,2);
  double xlmin = 0.45;
  double xlmax = 0.93;
  double ylmax = 0.88;
  double ylmin = ylmax - 0.07*(nprc+0.5);
  if ( ylmin < 0.40 ) ylmin = 0.40;
  string sprclab;
  //double ymax[4] = {0.002, 0.008, 0.15, 0.008};
  double ymax[4] = {0.0024, 0.0024, 0.0024, 0.0024};
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
      ph->GetYaxis()->SetTitle("Power/tick [(ke)^{2}/(20 kHz)]");
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
      ph->Scale(sfac);
      ph->SetLineColor(icol);
      ph->Print();
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
        poman->setTitle("DFT power");
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
      TLegend* pleg = poman->getLegend();
      pleg->SetMargin(0.10);  // Fraction used for symbol
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
        sspow << "#sqrt{#Sigma} = " << fixed << 1000.0*powrms << " e";
        slab += " " + sspow.str();
        pleg->AddEntry(pho, slab.c_str(), "l");
      }
    }
  }
  string sofnam = outpre + "dftmergelog_" + sprclab + "_run" + srun + ".{png,pdf}";
  cout << myname << "Printing " << sofnam << endl;
  omantop.print(sofnam);
  for ( int iman=0; iman<4; ++iman ) {
    TPadManipulator* pman = omantop.man(iman);
    pman->setLogY(false);
    pman->setMarginLeft(0.12);
    //TH1* phframe = pman->frameHist();
    //cout << "Frame X title: " << phframe->GetXaxis()->GetTitle() << endl;
  }
  sofnam = outpre + "dftmergelin_" + sprclab + "_run" + srun + ".{png,pdf}";
  omantop.update();
  cout << myname << "Printing " << sofnam << endl;
  omantop.print(sofnam);
}
