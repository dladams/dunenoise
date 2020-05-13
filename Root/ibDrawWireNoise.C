// ibDrawWireNoise.C
//
// David Adams
// May 2020
//
// Plot Iceberg wires with thickness indicating noise.
//
void ibDrawWireNoise(int run, string tfil, float nscal=10.0, float noff =0.8, float nmax =12.0) {
  vector<int> chwts;
  int ncha = 1280;
  if ( tfil.size() ) {
    chwts.resize(ncha);
    TPadManipulator* pman = TPadManipulator::read(tfil);
    if ( pman == nullptr ) {
      cout << "Unable to open " << tfil << endl;
      return;
    }
    TGraph* pg = pman->graph();
    double x, y;
    for ( int ipt=0; ipt<pg->GetN(); ++ipt ) {
      pg->GetPoint(ipt, x, y);
      int icha = x + 0.1;
      int wt = (y - noff)*nscal;
      if ( wt > nmax ) wt = nmax;
      if ( icha < ncha ) chwts[icha] = wt;
    }
  }
  string srun = to_string(run);
  for ( string view : {"z", "u", "v"} ) {
    for ( string side : {"front", "back"} ) {
      string fnam = "wirenoise/wirenoiseRun" + srun + view + side + ".png";
      string ttl = "Noise for run " + srun + " " + view + " " + side;
      ibPlotWires(view, side, fnam, ttl, chwts);
    }
  }
}
