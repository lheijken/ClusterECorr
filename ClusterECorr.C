#include <iostream>

#include <TQObject.h>
#include <RQ_OBJECT.h>
#include <TGClient.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>
#include <TLine.h>
#include <TF1.h>
#include <TGButton.h>
#include <TRootEmbeddedCanvas.h>
#include <TGNumberEntry.h>

using namespace std;

//-- The GUI class
class MyMainFrame : public TGMainFrame {
private:
  TRootEmbeddedCanvas *fEcanvas;
  TGNumberEntry       *fNumber;
  TH1D *h1D = 0;
  TH1D *hCorr = 0;
  TH1D *hCorrFinal = 0;

public:
  MyMainFrame(const TGWindow *p,UInt_t w,UInt_t h);
  virtual ~MyMainFrame();
  void NextSliceSkip();
  void NextSliceSave();
  void DrawSlice();
  void FitSliceLandau();
  void FitSliceGaus();
  void FitSliceCB();
  void Save();
  void SetFitRange();
};

//-- Global stuff
TH2D *h2D = 0;
int histslice = 0;
int slicewidth=2; // Multiples of existing x-bins
double fitrange_insigma = 1;
double crystalball_function_simple_highendtail(double x, double k, double sigma, double mean);
double crystalball_function_simple_highendtail(const double *x, const double *p);

//-- Class implementation
MyMainFrame::MyMainFrame(const TGWindow *p,UInt_t w,UInt_t h)
  : TGMainFrame(p,w,h) {

  //-- Create the histograms which the corrections will be stored in
  int nrcorrbins = (h2D->GetXaxis()->GetNbins())/slicewidth;
  double lowx = h2D->GetXaxis()->GetBinLowEdge(1);
  double upx = h2D->GetXaxis()->GetBinUpEdge(h2D->GetXaxis()->GetNbins());
  int extrabins = (h2D->GetXaxis()->GetNbins())%slicewidth;
  if(extrabins != 0){
    double newbinwidth = slicewidth*(h2D->GetXaxis()->GetBinWidth(1));
    nrcorrbins = nrcorrbins + extrabins;
    upx = newbinwidth*nrcorrbins;
  }
  hCorr = new TH1D("hCorrections","Extracted E_{k} corrections vs E_{k}^{rec}",nrcorrbins,lowx,upx); hCorr->Sumw2();
  double hstart = h2D->GetXaxis()->GetBinLowEdge(1);
  double hstop = h2D->GetXaxis()->GetBinUpEdge(h2D->GetXaxis()->GetNbins());
  int hbins = ((int)hstop-hstart) + 1;
  hCorrFinal = new TH1D("hCorrectionsFinal","Final E_{k} corrections vs E_{k}^{rec};E_{k}^{rec}",hbins,hstart,hstop);

  // Create canvas widget
  fEcanvas = new TRootEmbeddedCanvas("Ecanvas",this,w,h);
  AddFrame(fEcanvas, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 10,10,10,1));

  // Create a horizontal frame widget with buttons
  TGHorizontalFrame *hframe = new TGHorizontalFrame(this,w,50);

  TGTextButton *nextskip = new TGTextButton(hframe,"Skip (&n)");
  nextskip->Connect("Clicked()","MyMainFrame",this,"NextSliceSkip()");
  hframe->AddFrame(nextskip, new TGLayoutHints(kLHintsCenterX,5,5,3,4));

  TGTextButton *nextsave = new TGTextButton(hframe,"Save (&m)");
  nextsave->Connect("Clicked()","MyMainFrame",this,"NextSliceSave()");
  hframe->AddFrame(nextsave, new TGLayoutHints(kLHintsCenterX,5,5,3,4));

  TGTextButton *fitlandau = new TGTextButton(hframe,"FitLandau (&f)");
  fitlandau->Connect("Clicked()","MyMainFrame",this,"FitSliceLandau()");
  hframe->AddFrame(fitlandau, new TGLayoutHints(kLHintsCenterX,5,5,3,4));

  TGTextButton *fitgaus = new TGTextButton(hframe,"FitGaus (&g)");
  fitgaus->Connect("Clicked()","MyMainFrame",this,"FitSliceGaus()");
  hframe->AddFrame(fitgaus, new TGLayoutHints(kLHintsCenterX,5,5,3,4));

  TGTextButton *fitcb = new TGTextButton(hframe,"FitCB (&h)");
  fitcb->Connect("Clicked()","MyMainFrame",this,"FitSliceCB()");
  hframe->AddFrame(fitcb, new TGLayoutHints(kLHintsCenterX,5,5,3,4));

  fNumber = new TGNumberEntry(hframe, 0, 9,999, TGNumberFormat::kNESRealOne,TGNumberFormat::kNEANonNegative, TGNumberFormat::kNELLimitMinMax, 0, 999);
  fNumber->Connect("ValueSet(Long_t)", "MyMainFrame", this, "SetFitRange()");
  (fNumber->GetNumberEntry())->Connect("ReturnPressed()", "MyMainFrame", this, "SetFitRange()");
  hframe->AddFrame(fNumber, new TGLayoutHints(kLHintsCenterX,5,5,1,1));

  TGTextButton *save = new TGTextButton(hframe,"Save (&s)");
  save->Connect("Clicked()","MyMainFrame",this,"Save()");
  hframe->AddFrame(save, new TGLayoutHints(kLHintsCenterX,5,5,3,4));

  TGTextButton *exit = new TGTextButton(hframe,"Exit (&e)","gApplication->Terminate(0)");
  hframe->AddFrame(exit, new TGLayoutHints(kLHintsCenterX,5,5,3,4));

  AddFrame(hframe, new TGLayoutHints(kLHintsCenterX,2,2,2,2));

  // Set a name to the main frame
  SetWindowName("Simple Example");

  // Map all subwindows of main frame
  MapSubwindows();

  // Initialize the layout algorithm
  Resize(GetDefaultSize());

  // Map main frame
  MapWindow();
}

MyMainFrame::~MyMainFrame()
{
  // Clean up used widgets: frames, buttons, layout hints
  //Cleanup();
  delete this;
}

void MyMainFrame::NextSliceSkip()
{
  //-- Advance the slice index, check if we've reached the upper edge
  histslice++;
  if((histslice*slicewidth) > h2D->GetXaxis()->GetNbins())
    histslice--;

  //-- Delete the fit function
  TF1 *fitresult = 0;
  if(h1D) fitresult = h1D->GetFunction("fitfunc");
  if(fitresult){
    h1D->GetListOfFunctions()->Remove(fitresult);
  }
  //-- Find the bin range for the next slice
  int upbin = histslice*slicewidth;
  int lowbin = upbin - (slicewidth-1);
  double lowval = h2D->GetXaxis()->GetBinLowEdge(lowbin);  double upval = h2D->GetXaxis()->GetBinUpEdge(upbin);
  //-- Create the 1D distribution
  h1D = (TH1D*)h2D->ProjectionY(Form("Slice%d",histslice),lowbin,upbin);

  DrawSlice();
}

void MyMainFrame::NextSliceSave()
{
  //-- Save the fit results and then delete the fit function
  TF1 *fitresult = 0;
  if(h1D) fitresult = h1D->GetFunction("fitfunc");
  if(fitresult){
    //-- Find the current Ek value
    int upekbin = histslice*slicewidth;
    int lowekbin = upekbin - (slicewidth-1);
    double lowekval = h2D->GetXaxis()->GetBinLowEdge(lowekbin);
    double upekval = h2D->GetXaxis()->GetBinUpEdge(upekbin);
    double ekval = (lowekval + upekval)/2.;
    int ekbin = hCorr->GetXaxis()->FindFixBin(ekval);
    if(ekbin>(hCorr->GetXaxis()->GetNbins())) ekbin = hCorr->GetXaxis()->GetNbins();
    //-- Store the mean value and its error
    double meanres = fitresult->GetParameter(1);
    double meanerr = fitresult->GetParError(1);
    hCorr->SetBinContent(ekbin,meanres);
    hCorr->SetBinError(ekbin,meanerr);
    //-- Delete the fit function
    h1D->GetListOfFunctions()->Remove(fitresult);
    cout<<"Stored fit result for Ek = "<<ekval<<" MeV"<<endl;
  }
  else
    cout<<"No fit result to store"<<endl;


  //-- Advance the slice index, check if we've reached the upper edge
  histslice++;
  if((histslice*slicewidth) > h2D->GetXaxis()->GetNbins())
    histslice--;

  //-- Find the bin range for the next slice
  int upbin = histslice*slicewidth;
  int lowbin = upbin - (slicewidth-1);
  double lowval = h2D->GetXaxis()->GetBinLowEdge(lowbin);  double upval = h2D->GetXaxis()->GetBinUpEdge(upbin);
  //-- Create the 1D distribution
  h1D = (TH1D*)h2D->ProjectionY(Form("Slice%d",histslice),lowbin,upbin);

  DrawSlice();
}


void MyMainFrame::DrawSlice()
{
  //-- Find the range for the current slice
  int upbin = histslice*slicewidth;
  int lowbin = upbin - (slicewidth-1);
  double lowval = h2D->GetXaxis()->GetBinLowEdge(lowbin);
  double upval = h2D->GetXaxis()->GetBinUpEdge(upbin);
  //-- Create some lines to plot in the 2D distribution
  double ylow = h2D->GetYaxis()->GetBinLowEdge(1);
  double yup = h2D->GetYaxis()->GetBinUpEdge(h2D->GetYaxis()->GetNbins());
  TLine *llow = new TLine(lowval,ylow,lowval,yup);
  llow->SetLineWidth(2); llow->SetLineColor(2);
  TLine *lup = new TLine(upval,ylow,upval,yup);
  lup->SetLineWidth(2); lup->SetLineColor(2);
  //-- Get the fit function
  TF1 *fitresult = h1D->GetFunction("fitfunc");
  //-- Draw it all
  TCanvas *fCanvas = fEcanvas->GetCanvas();
  fCanvas->cd();
  if(histslice==1 && !fitresult) fCanvas->Divide(2,1);
  fCanvas->Modified();
  fCanvas->Update();
  fCanvas->cd(1);
  h2D->Draw("colz");
  llow->Draw(); lup->Draw();
  fCanvas->cd(2);
  h1D->Draw();
  if(fitresult) fitresult->DrawCopy("same");
  fCanvas->Update();
}

void MyMainFrame::FitSliceLandau()
{
  //-- Just in case they push the fit button before the Next one
  if(!h1D){
    cout<<"Go to the next slide first"<<endl;
    return;
  }

  double histlimlow = h1D->GetXaxis()->GetBinLowEdge(1);
  double histlimup =  h1D->GetXaxis()->GetBinUpEdge(h1D->GetXaxis()->GetNbins());

  //-- Fit function
  TF1 *f1Dfit = new TF1("fitfunc","[0]*TMath::Landau(x,[1],[2])",histlimlow,histlimup);
  f1Dfit->SetLineColor(2);

  //-- Get starting values from the 1D distribution
  double maxval = h1D->GetMaximum();
  double meanstart = h1D->GetXaxis()->GetBinCenter(h1D->GetMaximumBin());
  double widthstart = h1D->GetStdDev();
  //--- Adjust them to reasonable values
  if(maxval<0) maxval=1;
  if(meanstart<histlimlow || meanstart > histlimup) meanstart = (histlimup + histlimlow)/2.;
  if(widthstart > ((histlimup-histlimlow)/2.)) widthstart = ((histlimup-histlimlow)/2.);
  //-- Give them to the fit and impose limits
  f1Dfit->SetParameters(maxval,meanstart,widthstart);

  //-- Set the fit range, currently only based on a chosen multiple of the width around the mean value
  double fitlimlow = meanstart - (fitrange_insigma*widthstart);
  double fitlimup = meanstart + (fitrange_insigma*widthstart);
  f1Dfit->SetRange(fitlimlow,fitlimup);
  //-- Fit the 1D distribution
  h1D->Fit("fitfunc","B0R");

  DrawSlice();

}

void MyMainFrame::FitSliceGaus()
{
  //-- Just in case they push the fit button before the Next one
  if(!h1D){
    cout<<"Go to the next slide first"<<endl;
    return;
  }

  double histlimlow = h1D->GetXaxis()->GetBinLowEdge(1);
  double histlimup =  h1D->GetXaxis()->GetBinUpEdge(h1D->GetXaxis()->GetNbins());
  //-- Fit function
  TF1 *f1Dfit = new TF1("fitfunc","gaus",histlimlow,histlimup);
  f1Dfit->SetLineColor(2);

  //-- Get starting values from the 1D distribution
  double maxval = h1D->GetMaximum();
  double meanstart = h1D->GetXaxis()->GetBinCenter(h1D->GetMaximumBin());
  double widthstart = h1D->GetStdDev();
  //--- Adjust them to reasonable values
  if(maxval<0) maxval=1;
  if(meanstart<histlimlow || meanstart > histlimup) meanstart = (histlimup + histlimlow)/2.;
  if(widthstart > ((histlimup-histlimlow)/2.)) widthstart = ((histlimup-histlimlow)/2.);
  //-- Give them to the fit and impose limits
  f1Dfit->SetParameters(maxval,meanstart,widthstart);

  //-- Set the fit range, currently only based on a chosen multiple of the width around the mean value
  double fitlimlow = meanstart - (fitrange_insigma*widthstart);
  double fitlimup = meanstart + (fitrange_insigma*widthstart);
  f1Dfit->SetRange(fitlimlow,fitlimup);
  //-- Fit the 1D distribution
  h1D->Fit("fitfunc","B0R");

  DrawSlice();

}

void MyMainFrame::FitSliceCB()
{
  //-- Just in case they push the fit button before the Next one
  if(!h1D){
    cout<<"Go to the next slide first"<<endl;
    return;
  }

  double histlimlow = h1D->GetXaxis()->GetBinLowEdge(1);
  double histlimup =  h1D->GetXaxis()->GetBinUpEdge(h1D->GetXaxis()->GetNbins());
  //-- Fit function
  TF1 *f1Dfit = new TF1("fitfunc",crystalball_function_simple_highendtail,histlimlow,histlimup,4);
  f1Dfit->SetLineColor(2);

  //-- Get starting values from the 1D distribution
  double maxval = h1D->GetMaximum();
  double meanstart = h1D->GetXaxis()->GetBinCenter(h1D->GetMaximumBin());
  double widthstart = h1D->GetStdDev();
  //--- Adjust them to reasonable values
  if(maxval<0) maxval=1;
  if(meanstart<histlimlow || meanstart > histlimup) meanstart = (histlimup + histlimlow)/2.;
  if(widthstart > ((histlimup-histlimlow)/2.)) widthstart = ((histlimup-histlimlow)/2.);
  //-- Give them to the fit and impose limits
  f1Dfit->SetParameters(maxval,meanstart,widthstart,1);

  //-- Set the fit range, currently only based on a chosen multiple of the width around the mean value
  double fitlimlow = meanstart - (fitrange_insigma*widthstart);
  double fitlimup = meanstart + (fitrange_insigma*widthstart);
  f1Dfit->SetRange(fitlimlow,fitlimup);
  //-- Fit the 1D distribution
  h1D->Fit("fitfunc","B0R");

  DrawSlice();

}

void MyMainFrame::SetFitRange()
{
  //-- The fit range is (for now) set as multiple of the mean around the max value
  fitrange_insigma = fNumber->GetNumberEntry()->GetNumber();

}

void MyMainFrame::Save()
{
  //-- Create a smoothed version of the correction-histogram
  TH1D *hCorr_smooth = (TH1D*)hCorr->Clone("hCorr_smooth");
  //--- Search for the first non zero correction
  double firstcorrsm = 0;
  int ifirstcorrsm = 1;
  for(int bin=1; bin<=( hCorr->GetXaxis()->GetNbins() ); bin++){
    double currcorr = hCorr->GetBinContent(bin);
    if(currcorr != 0){
      firstcorrsm = currcorr;
      ifirstcorrsm = bin;
      break;
    }
  }
  //--- Then search for the last non zero correction
  double lastcorrsm = 0;
  int ilastcorrsm = hCorr_smooth->GetXaxis()->GetNbins();
  for(int bin=( hCorr_smooth->GetXaxis()->GetNbins() ); bin>0; bin--){
    double currcorr = hCorr_smooth->GetBinContent(bin);
    if(currcorr != 0){
      lastcorrsm = currcorr;
      ilastcorrsm = bin;
      break;
    }
  }
  //--- Then smooth the histogram
  int allbinssm = hCorr_smooth->GetXaxis()->GetNbins();
  hCorr_smooth->GetXaxis()->SetRange(ifirstcorrsm,ilastcorrsm);
  hCorr_smooth->Smooth(2,"R");
  hCorr_smooth->GetXaxis()->SetRange(1,allbinssm);

  //-- Fill the interpolation distribution with the smoothed corrected values
  for(int bin=1; bin<( hCorr_smooth->GetXaxis()->GetNbins() ); bin++){
    double currcorr = hCorr_smooth->GetBinContent(bin);
    if(currcorr != 0){
      int newbin = hCorrFinal->GetXaxis()->FindFixBin(hCorr_smooth->GetXaxis()->GetBinCenter(bin));
      hCorrFinal->SetBinContent(newbin,currcorr);
    }
  }
  //-- Fill out the empty bins in the correction histogram
  //   simple style: in the edges just use the edge value, inbetween interpolation
  //--- Search for the first non zero correction
  double firstcorr = 0;
  int ifirstcorr = 1;
  bool emptycorr = false;
  for(int bin=1; bin<=( hCorrFinal->GetXaxis()->GetNbins() ); bin++){
    double currcorr = hCorrFinal->GetBinContent(bin);
    if(currcorr != 0){
      firstcorr = currcorr;
      ifirstcorr = bin;
      break;
    }
    if(bin == hCorrFinal->GetXaxis()->GetNbins())
      emptycorr = true;
  }
  //--- Then search for the last non zero correction
  double lastcorr = 0;
  int ilastcorr = hCorrFinal->GetXaxis()->GetNbins();
  for(int bin=( hCorrFinal->GetXaxis()->GetNbins() ); bin>0; bin--){
    double currcorr = hCorrFinal->GetBinContent(bin);
    if(currcorr != 0){
      lastcorr = currcorr;
      ilastcorr = bin;
      break;
    }
  }
  //--- Then fill the adjusted correction histogram
  if(!emptycorr){
    double lowx = 0, upx = 0, lowy = 0, upy = 0, slope = 0, interc = 0;
    for(int bin=1; bin<=( hCorrFinal->GetXaxis()->GetNbins() ); bin++){
    //--- If we're below the first correction, fill with the first correction
      if(bin<ifirstcorr){
        hCorrFinal->SetBinContent(bin,firstcorr);
      }
      //--- If we're above the last correction, fill with the last correction
      else if(bin>ilastcorr){
        hCorrFinal->SetBinContent(bin,lastcorr);
      }
      //--- Otherwise, either just fill the correction or, if there is none, use an interpolation
      else{
        double currcorr = hCorrFinal->GetBinContent(bin);
        if(currcorr == 0){
          if(lowx==0){
            lowx = hCorrFinal->GetXaxis()->GetBinCenter(bin-1);
            lowy = hCorrFinal->GetBinContent(bin-1);
            int endbin = bin;
            while(hCorrFinal->GetBinContent(endbin) == 0){
              endbin++;
            }
            upx = hCorrFinal->GetXaxis()->GetBinCenter(endbin);
            upy = hCorrFinal->GetBinContent(endbin);
            slope = (upy-lowy)/(upx-lowx);
            interc = upy - slope*upx;
            currcorr = slope*(hCorrFinal->GetBinCenter(bin)) + interc;
            hCorrFinal->SetBinContent(bin,currcorr);
          }
          else{
            currcorr = slope*(hCorrFinal->GetBinCenter(bin)) + interc;
            hCorrFinal->SetBinContent(bin,currcorr);
          }
        }
        else{
          lowx = 0; upx = 0; lowy = 0; upy = 0; slope = 0; interc = 0;
        }
      }
    }
    int allbins = hCorrFinal->GetXaxis()->GetNbins();
    hCorrFinal->GetXaxis()->SetRange(ifirstcorr,ilastcorr);
    hCorrFinal->Smooth(2,"R");
    hCorrFinal->GetXaxis()->SetRange(1,allbins);
  }

  //-- Write the inputted and the extracted correction histograms to file
  TString inhname = h2D->GetName();
  TString outfname = "Corr_"+inhname+".root";
  TFile *fout = new TFile(outfname,"recreate");
  h2D->Write();
  hCorr->Write();
  hCorrFinal->Write();
  hCorr_smooth->Write();
  fout->Close();
}

double crystalball_function_simple_highendtail(double x, double k, double sigma, double mean)
{
  // evaluate the simple crystal ball function
  if (sigma < 0.)     return 0.;
  double z = (x - mean)/sigma;
  if (k < z){
    double e1 = (k*k)/2.;
    double e2 = k*z;
    return std::exp(e1-e2);
  }
  else {
    return std::exp(-0.5*z*z);
  }
}

double crystalball_function_simple_highendtail(const double *x, const double *p) {
  // if ((!x) || (!p)) return 0.; // just a precaution
  // [Constant] * ROOT::Math::crystalball_function(x, [k], [Sigma], [Mean])
  return (p[0] * crystalball_function_simple_highendtail(x[0], p[3], p[2], p[1]));
}

//-- Main function
void ClusterECorr()
{
  //-- Read in the file and extract the histogram
  TFile *fin = new TFile("INPUTFILE","read");
  if (fin->IsZombie()) {
    cout << "Error opening file" << endl;
    exit(-1);
  }
  else {
    h2D = (TH2D*)fin->Get("TrueRecCheck_ClusterE/h_Ek_TrueRecvsRec_em_CB");
    if(!h2D){
      cout << "Error opening histogram" << endl;
      exit(-1);
    }
  }

  //-- Start the gui
  new MyMainFrame(gClient->GetRoot(),800,400);

}
