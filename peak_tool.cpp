#include "peak_tool.hpp"

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
 {
    string peakFileName="";
    if(argc>1){
      peakFileName=argv[1];
    }
    else{
      cerr<<endl<<"Input should be a BED file with ChIP peaks (preferably peak-summit location)" \
	  <<endl<<"Example of 1st 5 columns in BED file " \
	  <<endl<<"chr1    762714  762715  MACS_peak_7     439.21";
      peakFileName="test.bed";
      WAITUSER
    }
    ifstream::pos_type fileSize=get_file_size(peakFileName);
    cerr<<endl<<"Using input file:"<<peakFileName;
    cerr<<endl<<"BED file size "<<fileSize<<" bytes";
    ifstream peakFile;
    peakFile.open(peakFileName.c_str());
    string line;
    vector <string> lines;
    while(getline(peakFile,line)){
      lines.push_back(line);
    }
    for(vector<string>::iterator l=lines.begin();l!=lines.end();++l){
      string bedLine=*l;
      if(bedLine.find_first_of("#")!=string::npos)//remove commented lines
	continue;
      stringstream ssLine(*l);
      vector<string> bedElements; 
      unsigned countColumns=0;
      while(ssLine.good()){
	  string element=".";
	  ssLine>>element;
	  ++countColumns;
	  //cout<<endl<<element<<"\t"<<countColumns;
	  bedElements.push_back(element);
      }
      if(bedElements.size()>=3){
	string s_chrNum=bedElements[0];
	string s_chrStart=bedElements[1];
	string s_chrEnd=bedElements[2];
	string name=".";
	string s_score=".";
	if(bedElements.size()>=5){
	  name=bedElements[3];
	  s_score=bedElements[4];
	}
	unsigned chrNum=get_chr_number(s_chrNum);
	if(chrNum==0){
	  continue;
	}
	unsigned chrStart=0;
	stringstream ss_chrStart(s_chrStart);
	ss_chrStart>>chrStart;
	unsigned chrEnd=0;
	stringstream ss_chrEnd(s_chrEnd);
	ss_chrEnd>>chrEnd;
	float score=0.0;
	stringstream ss_score(s_score);
	ss_score>>score;

	unsigned peakMid=chrStart/2+chrEnd/2;
	vector<unsigned short int> promoterIndexes;
	promoterIndexes.push_back(promoterMapsW[chrNum].featureMapChromosome[peakMid]);
	promoterIndexes.push_back(promoterMapsC[chrNum].featureMapChromosome[peakMid]);
// 	unsigned short int promoterIndexTestW=promoterMapsW[chrNum].featureMapChromosome[peakMid];
// 	unsigned short int promoterIndexTestC=promoterMapsC[chrNum].featureMapChromosome[peakMid];
	unsigned short int transcriptIndexTest=transcriptMaps[chrNum].featureMapChromosome[peakMid];
	bool enhancerTest=(enhancerBitmaps[chrNum].featureBitmapChromosome[peakMid]);
	bool superenhancerTest=(superenhancerBitmaps[chrNum].featureBitmapChromosome[peakMid]);
	bool promoterFound=0;
	for(vector<unsigned short int>::iterator i=promoterIndexes.begin(); i!=promoterIndexes.end(); ++i){
	  if((*i)>0){
	    promoterFound=1;
	    cout<<endl<<"chr"<<get_chr_text(chrNum)<<"\t"<<peakMid<<"\t"<<name<<"\t"<<score;
	    cout<<"\tPROMOTER";
	    output_gene_data(chrNum, (*i), peakMid, transcriptsData);
	  }
	}
	if (promoterFound==0){
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
	    output_gene_data(chrNum, transcriptIndexTest, peakMid, transcriptsData);
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
	    cout<<"\t.\t.\t.\t.\t.\t.";
	  }
	}
      }
    }
    lines.clear();
  }
  exonsData.clear();
  genesData.clear();
  transcriptsData.clear();
  exonBitmaps.clear();
  enhancerBitmaps.clear();
  promoterMapsW.clear();
  promoterMapsC.clear();
  transcriptMaps.clear();
  return 0;
  
}