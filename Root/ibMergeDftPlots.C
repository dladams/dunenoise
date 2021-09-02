// ibMergeDftPlots.C
//
// David Adams
// May 2020
//
// Merge Iceberg before and after CNR plots.

TPadManipulator* ibMergeDftPlots(string infilPat, string srun, string outfilPat, string splots) {
  const string myname = "ibMergeDftPlots: ";
  using Index = unsigned int;
  LineColors lc;
  string sfrun = srun;
  string sfrun1;
  string sfrun2;
  string::size_type ipos = sfrun.find("-");
  if ( ipos == string::npos ) {
    while ( sfrun.size() < 6 ) sfrun = "0" + sfrun;
  } else {
    sfrun1 = sfrun.substr(0, ipos);
    sfrun2 = sfrun.substr(ipos+1);
    while ( sfrun1.size() < 6 ) sfrun1 = "0" + sfrun1;
    while ( sfrun2.size() < 6 ) sfrun2 = "0" + sfrun2;
    sfrun = sfrun1 + "-" + sfrun2;
  }
  vector<string> srecs = {"tai", "cnr"};     // Iceberg 4
  vector<string> labs = {"Before CNR", "After CNR", ""};
  vector<string> srecs5cal = {"cal", "cnr"};
  vector<string> labs5cal = {"Calibrated", "After CNR", ""};
  vector<string> srecs5ped = {"ped", "pnr"};
  vector<string> labs5ped = {"After PUP", "After CNR", ""};
  if ( splots == "cal:cnr" ) {
    srecs = srecs5cal;
    labs = labs5cal;
  } else if ( splots == "ped:pnr" ) {
    srecs = srecs5ped;
    labs = labs5ped;
  } else if ( splots != "tai:cnr" ) {
    cout << myname << "Invalid splots: " << splots << endl;
    return nullptr;
  }
  vector<int> icols = {lc.blue(), lc.red()};
  vector<int> lwids = {3, 2};
  vector<float> ylabs = {0.73, 0.66};
  vector<string> sttls = {"#bf{DUNE:Iceberg}", "", "Run " + srun};
  vector<float> yttls = {0.40, 0.34, 0.28};
  vector<string> planes = {"Z1", "U", "Z2", "V"};
  TPadManipulator* pmanout = nullptr;
  Index nman = 4;
  pmanout = new TPadManipulator(1400, 1000);
  pmanout->split(2, 2);
  float tsiz = 0.05;
  string sttl0 = "DFT power for run " + srun + " plane ";
  float xlab = 0.50;
  float xttl = 0.60;
  for ( Index irec=0; irec<srecs.size(); ++irec ) {
    string srec = srecs[irec];
    int icol = icols[irec];
    StringManipulator sman(infilPat, true);
    sman.replace("%RECO%", srec);
    sman.replace("%RUN%", sfrun);
    sman.replace("%RUN1%", sfrun1);
    sman.replace("%RUN2%", sfrun2);
    string sfilin = sman.str();
    TPadManipulator* pman = TPadManipulator::read(sfilin);
    if ( pman == nullptr ) {
      cout << myname << "Unable to open " << sfilin << endl;
      //return nullptr;
      continue;
    }
    string hopt = irec ? "hist same" : "hist";
    for ( Index iman=0; iman<nman; ++iman ) {
      cout << "====== " << irec << "-" << iman << endl;
      TPadManipulator* pmano = pmanout->man(iman);
      TPadManipulator* pmani = pman->man(iman);
      TH1* ph = pmani->hist();
      ph->Scale(1.e6);
      string sttly = "Power/tick/channel [e^{2}]";
      ph->GetYaxis()->SetTitle(sttly.c_str());
      ph->SetLineColor(icol);
      ph->SetLineWidth(lwids[irec]);
      ph->SetTitle("");
      // Extract the noise in this view and the total.
      float hsum = ph->Integral();
      float hsumo = ph->Integral(0, ph->GetNbinsX()+1);
      float nsum = hsum > 0.0 ? sqrt(hsum) : 0.0;
      float nsumo = hsumo > 0.0 ? sqrt(hsumo) : 0.0;
      string snsum = to_string(int(nsum + 0.5));
      string snsumo = to_string(int(nsumo + 0.5));
      pmano->add(ph, hopt);
      vector<TLatex*> txts;
      for ( TObject* pobj : pmani->objects() ) {
        TLatex* ptxt = dynamic_cast<TLatex*>(pobj);
        cout << pobj->GetTitle() << endl;
        if ( ptxt != nullptr ) {
          txts.push_back(ptxt);
        }
      }
      string slab = labs[irec] + " ";
      slab += "#sqrt{#Sigma} = " + snsum + " e";
      TLatex* ptxt = new TLatex(xlab, ylabs[irec], slab.c_str());
      ptxt->SetNDC();
      ptxt->SetTextColor(icol);
      ptxt->SetTextFont(42);
      ptxt->SetTextSize(tsiz);
      pmano->add(ptxt, "");
    }
  }
  for ( Index iman=0; iman<nman; ++iman ) {
    TPadManipulator* pman = pmanout->man(iman);
    for ( Index ittl=0; ittl<sttls.size(); ++ittl ) {
      string sttl = sttls[ittl];
      if ( sttl.size() == 0 ) sttl = planes[iman] + " planes";;
      TLatex* ptxt = new TLatex(xttl, yttls[ittl], sttl.c_str());
      ptxt->SetNDC();
      ptxt->SetTextFont(42);
      ptxt->SetTextSize(tsiz);
      pman->add(ptxt, "");
    }
    pman->setLabelSizeX(tsiz);
    pman->setLabelSizeY(tsiz);
    pman->setTitleSize(tsiz);
    if ( false ) {
      pman->setMarginLeft(0.13);
      pman->setMarginRight(0.04);
    } else {
      pman->setMarginLeft(0.14);
      pman->setMarginTop(0.05);
    }
    pman->setMarginBottom(0.12);
    pman->setRangeY(0, 300.0);
    pman->setTitle("");
    pman->addAxis();
    pman->centerAxisTitles();
  }
  StringManipulator sman(outfilPat, true);
  sman.replace("%RUN%", srun);
  sman.replace(":", "_");
  string sfilout = sman.str();
  cout << myname << "Printing " << sfilout << endl;
  pmanout->print(sfilout);
  return pmanout;
}

