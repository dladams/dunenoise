void pdMergeDftPlots(string fins, string sonam, string slegin="After tail removal:,After CNR:", string slabin ="", double ylogmax =0.2) {
  string myname = "pdMergeDftPlots: ";
  StringManipulator smfins(fins, false);
  StringManipulator::StringVector finv = smfins.split(",");
  int nfin = finv.size();
  cout << myname << "File count is " << nfin << endl;
  StringManipulator smlegs(slegin, false);
  StringManipulator::StringVector slegs = smlegs.split(",");
  int nleg = slegs.size();
  cout << myname << "Label count is " << nleg << endl;
  TPadManipulator omantop(1400, 1000);
  omantop.split(2,2);
  double xlmin = 0.40;
  double xlmax = 0.93;
  double ylmax = 0.90;
  double ylmin = ylmax - 0.06*(nfin+0.5);
  if ( ylmin < 0.40 ) ylmin = 0.40;
  //double ymax[4] = {0.002, 0.008, 0.15, 0.008};
  double ymax[4] = {0.0033, 0.009, 0.18, 0.009};
  vector<string> padlabs = {"C wires", "U wires", "Z wires", "V wires"};
  for ( int ifin=0; ifin<nfin; ++ifin ) {
    string fnam = finv[ifin];
    string dopt = "HIST";
    int icol = LineColors::color(ifin);
    if ( ifin ) dopt += " SAME";
    TPadManipulator* pmantop = TPadManipulator::read(fnam);
    if ( pmantop == nullptr ) {
      cout << myname << "Unable to find tpm in " << fnam << endl;
      return;
    }
    for ( int iman=0; iman<4; ++iman ) {
      TPadManipulator* pman = pmantop->man(iman);
      TH1* ph = (TH1*) pman->hist()->Clone();
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
      if ( ifin == 0 ) {
        poman->addLegend(xlmin, ylmin, xlmax, ylmax);
      }
      cout << myname << "Extra object count for pad " << iman << ": " << poman->objects().size() << endl;
      if ( ifin + 1 == nfin ) {
        poman->setLogY();
        poman->setRangeY(0.0, ymax[iman]);
        poman->setLogRangeY(1.e-5, ylogmax);
        poman->setRangeY(0, ylogmax);
        poman->addAxis();
        //poman->setGrid();
        poman->setTitle("DFT power");
        double xlab = 0.65;
        double ylab = slabin.size() ? 0.28 : 0.22;
        TLatex* ptxt = new TLatex(xlab, ylab, "#bf{ProtoDUNE-SP}");
        ptxt->SetNDC();
        ptxt->SetTextFont(42);
        poman->add(ptxt);
        if ( slabin.size() ) {
          ylab -= 0.06;
          ptxt = new TLatex(xlab, ylab, slabin.c_str());
          ptxt->SetNDC();
          ptxt->SetTextFont(42);
          poman->add(ptxt);
        }
        ylab -= 0.06;
        ptxt = new TLatex(xlab, ylab, padlabs[iman].c_str());
        ptxt->SetNDC();
        ptxt->SetTextFont(42);
        poman->add(ptxt);
      }
      poman->getYaxis()->SetTitle("Power [(ke)^{2}/channel/tick]");
      TLegend* pleg = poman->getLegend();
      pleg->SetMargin(0.15);  // Fraction used for symbol
      if ( pleg == nullptr ) {
        cout << myname << "ERROR: Legend not found." << endl;
      } else {
        string slab = ifin < nleg ? slegs[ifin] : "";
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
  string sofnam = "outfiles/dftmergelog" + sonam + ".{png,pdf,tpad}";
  cout << myname << "Printing " << sofnam << endl;
  omantop.print(sofnam);
  for ( int iman=0; iman<4; ++iman ) {
    TPadManipulator* pman = omantop.man(iman);
    pman->setLogY(false);
    pman->setMarginLeft(0.12);
    //TH1* phframe = pman->frameHist();
    //cout << "Frame X title: " << phframe->GetXaxis()->GetTitle() << endl;
  }
  sofnam = "outfiles/dftmergelin" + sonam + ".{png,pdf,tpad}";
  omantop.update();
  cout << myname << "Printing " << sofnam << endl;
  omantop.print(sofnam);
}
