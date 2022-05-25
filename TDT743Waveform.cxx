#include "TDT743Waveform.h"

#include "TDT743RawData.hxx"
#include "TDirectory.h"


/// Reset the histograms for this canvas
TDT743Waveform::TDT743Waveform(){

  SetSubTabName("DT743 Waveforms");
  SetUpdateOnlyWhenPlotted(true);
  
  CreateHistograms();
  FrequencySetting = -1;
}


void TDT743Waveform::CreateHistograms(){

  // check if we already have histogramss.
  char tname[100];
  sprintf(tname,"DT743_%i",0);

  TH1D *tmp = (TH1D*)gDirectory->Get(tname);
  if (tmp) return;

  int fWFLength = 7; // Need a better way of detecting this...
  int numSamples = fWFLength / 1;
  
  // Otherwise make histograms
  clear();
  
  for(int i = 0; i < 8; i++){ // loop over 2 channels
    
    char name[100];
    char title[100];
    sprintf(name,"DT743_%i",i);
    
    sprintf(title,"DT743 Waveform for channel=%i",i);	
    
    TH1D *tmp = new TH1D(name, title, numSamples, 0, fWFLength);
    tmp->SetXTitle("ns");
    tmp->SetYTitle("ADC value");
    
    push_back(tmp);
  }
  std::cout << "TDT743Waveform done init...... " << std::endl;
  FrequencySetting = -1;
}


void TDT743Waveform::UpdateHistograms(TDataContainer& dataContainer){
  
  int eventid = dataContainer.GetMidasData().GetEventId();
  int timestamp = dataContainer.GetMidasData().GetTimeStamp();
  
  TDT743RawData *dt743 = dataContainer.GetEventData<TDT743RawData>("43FS");
  
  if(dt743){      
    
    
    std::vector<RawChannelMeasurement> measurements = dt743->GetMeasurements();

    bool changeHistogram = false; 
    for(int i = 0; i < measurements.size(); i++){
           
      int chan = measurements[i].GetChannel();
      int nsamples = measurements[i].GetNSamples();

      // Check if we need to change timebase and number of samples using first channel.
      if(i == 0){
	if(FrequencySetting == -1 or FrequencySetting != measurements[i].GetFrequency()
	   or nsamples != GetHistogram(chan)->GetNbinsX()){
	  std::cout << "New setting for histogram : " << FrequencySetting << " " << measurements[i].GetFrequency()
		    << " " << nsamples << " " << GetHistogram(chan)->GetNbinsX() << std::endl;
	  changeHistogram = true;
	  FrequencySetting = measurements[i].GetFrequency();
	}  	  
      }
      
      // Update the histogram
      if(changeHistogram){
	float ns_per_bin = 0;
	if(FrequencySetting == 0) ns_per_bin = 0.3125;
	else if(FrequencySetting == 1) ns_per_bin = 0.625;
	else if(FrequencySetting == 2) ns_per_bin = 1.25;
	else if(FrequencySetting == 3) ns_per_bin = 2.5;
	GetHistogram(chan)->Reset();
	GetHistogram(chan)->SetBins(nsamples, 0, nsamples*ns_per_bin);
      }
      
      //std::cout << "Nsamples " <<  measurements[i].GetNSamples() << std::endl;
      for(int ib = 0; ib < nsamples; ib++){
	GetHistogram(chan)->SetBinContent(ib+1, measurements[i].GetSample(ib));
      }
      
    }
  }
  
}



void TDT743Waveform::Reset(){
  
  
  for(int i = 0; i < 8; i++){ // loop over channels
    int index =  i;
    
    // Reset the histogram...
    for(int ib = 0; ib < GetHistogram(index)->GetNbinsX(); ib++) {
      GetHistogram(index)->SetBinContent(ib, 0);
    }
    
    GetHistogram(index)->Reset();
    
  }
}
