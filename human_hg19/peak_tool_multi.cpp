#include "peak_tool_multi.hpp"

int main (int argc, char *argv[]){  
  
  vector<feature_bitmaps> enhancerBitmaps;
  string enhancerFileName="/home/amitra/proj/abhishekNGS/huge_files/human/enhancers/enhancers.bed";
  int ret;
  ret=load_enhancers(enhancerBitmaps,enhancerFileName);
  if (ret==0){
    enhancerFileName="enhancers.bed";
    ret=load_enhancers(enhancerBitmaps,enhancerFileName);
  }
 
  vector<feature_bitmaps> superenhancerBitmaps;
  string superenhancerFileName="/home/amitra/proj/abhishekNGS/huge_files/human/enhancers/super-enhancers.bed";
  ret=load_enhancers(superenhancerBitmaps,superenhancerFileName);
  if (ret==0){
    enhancerFileName="enhancers.bed";
    ret=load_enhancers(enhancerBitmaps,enhancerFileName);
  }

  if(argc<2){
    cerr << endl << "peaks_tool_multi <peak file primary> <peak file 2>...<peak file n>"<<endl;;
    return 0;
  }
  vector<string> bedFileNames;
  vector <bedFilePeak> peakFiles;
  for(int nFiles=1;nFiles<argc;++nFiles){
      cerr<<endl<<"File"<<nFiles<<"\t"<<argv[nFiles];
      bedFileNames.push_back(argv[nFiles]);
  }
  cerr<<endl<<"Primary Peak File = "<<bedFileNames[0];
  unsigned numBedFiles=bedFileNames.size();
  cerr<<endl<<"Total files = "<<numBedFiles;
  read_peak_files(bedFileNames, peakFiles);
  cerr<<endl<<"Finding Peak Neighbors";
  primaryPeakData primaryPeak;
  find_primary_peak_neighbors(primaryPeak, peakFiles, numBedFiles);
 
  cerr<<endl<<"Reading Gencode";
  vector<feature_bitmaps> exonBitmaps;
  vector<feature_maps> promoterMapsW;
  vector<feature_maps> promoterMapsC;
  vector<feature_maps> transcriptMaps;
  vector<genes> genesData;
  vector<transcripts> transcriptsData;
  vector<exons> exonsData;
  string gencodeFileName="/home/amitra/proj/abhishekNGS/huge_files/human/gencode/gencode.v19.annotation.gtf";
  ret=read_gencode_annotation_file(exonBitmaps, promoterMapsW, promoterMapsC, transcriptMaps, \
				genesData, transcriptsData, exonsData, gencodeFileName);
  if(ret==0){
    gencodeFileName="gencode.v19.annotation.gtf";
    ret=read_gencode_annotation_file(exonBitmaps, promoterMapsW, promoterMapsC, transcriptMaps, \
				genesData, transcriptsData, exonsData, gencodeFileName);
  }
  if(ret==0){
    cerr<<endl<<"Cannot annotate, because Gencode file is not loaded.";
    return 0;
  }
    
  cerr<<endl<<"Annotating";
  unsigned int peakCount=0;
  for(vector<bedData5>::iterator p=primaryPeak.peakList.begin(); p!=primaryPeak.peakList.end(); ++p,++peakCount){
      unsigned int peakMid=(*p).chrStart/2 + (*p).chrEnd/2;
      unsigned int chrNum=(*p).chrNum;
      string name=(*p).name;
      float score=(*p).score;
      vector<bedData5> neighborPeakList=primaryPeak.neighborData[peakCount].neighborPeakList;
      vector<unsigned short int> promoterIndexes;
      promoterIndexes.push_back(promoterMapsW[chrNum].featureMapChromosome[peakMid]);
      promoterIndexes.push_back(promoterMapsC[chrNum].featureMapChromosome[peakMid]);
//       unsigned short int promoterIndexTestW=promoterMapsW[chrNum].featureMapChromosome[peakMid];
//       unsigned short int promoterIndexTestC=promoterMapsC[chrNum].featureMapChromosome[peakMid];
      unsigned short int transcriptIndexTest=transcriptMaps[chrNum].featureMapChromosome[peakMid];
      bool enhancerTest=enhancerBitmaps[chrNum].featureBitmapChromosome[peakMid];
      bool superenhancerTest=superenhancerBitmaps[chrNum].featureBitmapChromosome[peakMid];
      bool promoterFound=0;
      if(primaryPeak.numNeighbors[peakCount]<numBedFiles-1)continue;//only output primary peak with all neighbors
      for(vector<unsigned short int>::iterator i=promoterIndexes.begin(); i!=promoterIndexes.end(); ++i){
	if((*i)>0){
	  promoterFound=1;
	  cout<<endl<<"chr"<<get_chr_text(chrNum)<<"\t"<<peakMid<<"\t"<<name<<"\t"<<score;
	  cout<<"\tPROMOTER";
	  output_gene_data(chrNum, (*i), peakMid, transcriptsData, neighborPeakList);
	}
      }
      if(promoterFound==0){
	if(transcriptIndexTest>0){
	    cout<<endl<<"chr"<<get_chr_text(chrNum)<<"\t"<<peakMid<<"\t"<<name<<"\t"<<score;
	    bool exonTest=0;
	    exonTest=exonBitmaps[chrNum].featureBitmapChromosome[peakMid];
	    if(exonTest){
	      cout<<"\tEXON";
	    }
	    else{
	      cout<<"\tINTRON";
	    }
    	  output_gene_data(chrNum, transcriptIndexTest, peakMid, transcriptsData, neighborPeakList);
	}
	else{
	  cout<<endl<<"chr"<<get_chr_text(chrNum)<<"\t"<<peakMid<<"\t"<<name<<"\t"<<score;
	  if(enhancerTest>0||superenhancerTest>0){
	    cout<<"\tENHANCER";	
	    if(superenhancerTest>0)cout<<"-SUPER";
	  }
	  else{
	   cout<<"\tINTERGENIC";
	  }
	  cout<<"\t.\t.\t.\t.\t.\t.";//no transcript or promoter data
	  for(vector<bedData5>::iterator n=neighborPeakList.begin(); n!=neighborPeakList.end();++n){
	      unsigned peakMid2=(*n).chrEnd/2+(*n).chrStart/2;
	      if(peakMid2!=0){
		cout<<"\t"<<peakMid2<<"\t"<<(*n).name<<"\t"<<(*n).score<<"\t"<<(int)peakMid2-(int)peakMid<<"\t.";
	      }
	      else{
		cout<<"\t.\t.\t.\t.\t.";
	      }
	  }
	}
      }
    }
  exonsData.clear();
  genesData.clear();
  transcriptsData.clear();
  exonBitmaps.clear();
  enhancerBitmaps.clear();
  promoterMapsW.clear();
  promoterMapsC.clear();
  transcriptMaps.clear();
  bedFileNames.clear();
  peakFiles.clear();
  return 0;
}