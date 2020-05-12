// ibMergeDftPlots.C
//
// David Adams
// May 2020
//
// Merge Iceberg before and after CNR plots.

TPadManipulator* ibMergeDftPlots(string infilPat, int irun, string outfilPat) {
  const string myname = "ibMergeDftPlots: ";
  using Index = unsigned int;
  LineColors lc;
  string srun = to_string(irun);
  string sfrun = srun;
  while ( sfrun.size() < 6 ) sfrun = "0" + sfrun;
  vector<string> srecs = {"tai", "cnr"};
  vector<int> icols = {lc.blue(), lc.red()};
  vector<int> lwids = {3, 2};
  vector<string> labs = {"Before CNR", "After CNR"};
  vector<float> ylabs = {0.80, 0.73};
  vector<string> planes = {"z1", "u", "z2", "v"};
  TPadManipulator* pmanout = nullptr;
  Index nman = 4;
  pmanout = new TPadManipulator(1400, 1000);
  pmanout->split(2, 2);
  float tsiz = 0.045;
  string sttl0 = "DFT power for run " + srun + " plane ";
  for ( Index irec=0; irec<srecs.size(); ++irec ) {
    string srec = srecs[irec];
    int icol = icols[irec];
    StringManipulator sman(infilPat, true);
    sman.replace("%RECO%", srec);
    sman.replace("%RUN%", sfrun);
    string sfilin = sman.str();
    TPadManipulator* pman = TPadManipulator::read(sfilin);
    if ( pman == nullptr ) {
      cout << myname << "Unable to open " << sfilin << endl;
      return nullptr;
    }
    string hopt = irec ? "hist same" : "hist";
    float xlab = 0.50;
    for ( Index iman=0; iman<nman; ++iman ) {
      cout << "====== " << irec << "-" << iman << endl;
      TPadManipulator* pmano = pmanout->man(iman);
      TPadManipulator* pmani = pman->man(iman);
      TH1* ph = pmani->hist();
      ph->Scale(1.e6);
      string sttly = ph->GetYaxis()->GetTitle();
      sttly += " [e^{2}]";
      ph->GetYaxis()->SetTitle(sttly.c_str());
      ph->SetLineColor(icol);
      ph->SetLineWidth(lwids[irec]);
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
      string sttl = sttl0 + planes[iman];
      pmano->setTitle(sttl);
    }
  }
  for ( Index iman=0; iman<nman; ++iman ) {
    TPadManipulator* pman = pmanout->man(iman);
    pman->setLabelSizeX(tsiz);
    pman->setLabelSizeY(tsiz);
    pman->setTitleSize(tsiz);
    pman->setMarginLeft(0.13);
    pman->setMarginRight(0.04);
    pman->setMarginBottom(0.12);
    pman->setRangeY(0, 300.0);
    pman->addAxis();
  }
  StringManipulator sman(outfilPat, true);
  sman.replace("%RUN%", srun);
  string sfilout = sman.str();
  cout << myname << "Printing " << sfilout << endl;
  pmanout->print(sfilout);
  return pmanout;
}

