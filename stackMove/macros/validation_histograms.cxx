// ROOT macro producing histograms from MCStepAnalysis

#include <iostream>
#include <vector>
#include <string>

#include "TH1D.h"
#include "TCanvas.h"
#include "TLegend.h"

#include "MCStepLogger/MCAnalysisFileWrapper.h"
#include "MCStepLogger/MetaInfo.h"


const std::vector<int> colours = {kRed+2, kBlue+3, kGreen+3, kYellow+3, kCyan+3};
const std::vector<int> lineStyles = {1, 7, 3};


// Produce overlay histograms for some observables and pass to canvases
void validationHistograms(const std::vector<std::string>& inFilePaths,
                          const std::vector<std::string>& observables,
                          const std::vector<std::string>& xAxisLabels,
                          const std::vector<std::string>& yAxisLabels,
                          const std::string& outDir)
{


  if(observables.size() != xAxisLabels.size() ||
     xAxisLabels.size() != yAxisLabels.size()) {
    std::cerr << "Number of observabled and axis labels must have "
              << "the same size.\n";
    return;
  }
  std::vector<o2::mcstepanalysis::MCAnalysisFileWrapper> fileWrappers;
  // Need as amny wrappers as file paths.
  fileWrappers.reserve(inFilePaths.size());

  // Open all files in an MCAnalysisFileWrapper
  for(int i = 0; i < inFilePaths.size(); i++) {
    o2::mcstepanalysis::MCAnalysisFileWrapper wrapper;
    wrapper.read(inFilePaths[i]);
    fileWrappers.push_back(wrapper);
  }

  // Helper vector to make sure the histograms are properly visible
  std::vector<float> maxima(observables.size(), 0.);
  for(int i = 0; i < observables.size(); i++) {
    for(int j = 0; j < fileWrappers.size(); j++) {
      // Note that this will fail in case the histogram is not there.
      TH1D& histo = fileWrappers[j].getHistogram<TH1D>(observables[i]);
      // update the maximum
      if(histo.GetMaximum() > maxima[i]) {
        maxima[i] = histo.GetMaximum();
      }
    }
  }
  TCanvas* canvas = nullptr;
  TLegend* legend = nullptr;
  // Now fill the histograms by looping over file wrappers and observables
  for(int i = 0; i < observables.size(); i++) {
    // Create new canvas and legend for a new observable.
    legend = new TLegend(0.7, 0.7, 0.9, 0.9);
    canvas = new TCanvas(observables[i].c_str(), "", 600, 500);
    std::string outFilePath(outDir + "/" + observables[i]);
    canvas->cd();
    for(int j = 0; j < fileWrappers.size(); j++) {
      // Note that this will fail in case the histogram is not there.
      TH1D& histo = fileWrappers[j].getHistogram<TH1D>(observables[i]);
      histo.SetMaximum(1.2 * maxima[i]);
      histo.GetXaxis()->SetTitle(xAxisLabels[i].c_str());
      histo.GetYaxis()->SetTitle(yAxisLabels[i].c_str());
      histo.SetLineColor(colours[j%colours.size()]);
      histo.SetLineStyle(lineStyles[j%lineStyles.size()]);
      histo.SetLineWidth(2);
      histo.SetStats(kFALSE);
      histo.Draw("same hist");
      legend->AddEntry(&histo, fileWrappers[j].getAnalysisMetaInfo().label.c_str());
    }
    legend->Draw("same");
    canvas->SaveAs(std::string(outFilePath + ".eps").c_str());
    canvas->SaveAs(std::string(outFilePath + ".png").c_str());
  }
  std::cerr << "All histograms produced\n";
}

int main(int argc, char* argv[])
{
  if(argc < 3) {
    std::cerr << "The first argument is assumed to be the output directory "
              << "(which msut exist already). "
              << "Furthermore, at least one analysis file is required\n";
    return 1;
  }
  const std::string outDir(argv[1]);

  std::vector<std::string> inFilePaths;
  for(int i = 2; i < argc; i++) {
    inFilePaths.push_back(argv[i]);
  }
  // Histogram names
  std::vector<std::string> observables = {"stepsXPerEvent",
                                          "stepsYPerEvent",
                                          "stepsZPerEvent",
                                          "nStepsPerEvent",
                                          "meanStepSizePerVolPerEvent",
                                          "relNStepsPerPDGPerEvent",
                                          "relNStepsPerVolPerEvent",
                                          "nSecondariesPerVolPerEvent",
                                          "stepsEnergyPerEvent",
                                          "nSteps"};
  // x axis labels
  std::vector<std::string> xAxisLabels = {"x [cm]",
                                          "y [cm]",
                                          "z [cm]",
                                          "",
                                          "volume",
                                          "PDG",
                                          "volume",
                                          "volume",
                                          "energy [GeV]",
                                          ""};
  // y axis labels
  std::vector<std::string> yAxisLabels = {"#counts/event",
                                          "#counts/event",
                                          "#counts/event",
                                          "#steps/event",
                                          "step size [cm]",
                                          "1/event",
                                          "1/event",
                                          "#secondaries/event",
                                          "#counts/event",
                                          "#steps"};


  validationHistograms(inFilePaths, observables, xAxisLabels, yAxisLabels, outDir);
  return 0;
}
