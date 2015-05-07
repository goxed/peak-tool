#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#define WAITUSER    cerr<<endl<<"press return"; cin.ignore();
#define X_CHROMOSOME 1000
#define Y_CHROMOSOME 1001
#define M_CHROMOSOME 1002

#define PROMOTER_UPSTREAM 2000
#define PROMOTER_DOWNSTREAM 500

using namespace std;

struct transcript
{
  unsigned chrNum;
  string geneName;
  string transcriptName;
  unsigned chrStart;
  unsigned chrEnd;
  bool WC;
  string transcriptType;
};

struct exon
{
  unsigned chrNum;
  string geneName;
  unsigned exonNumber;
  string transcriptName;
  unsigned chrStart;
  unsigned chrEnd;
  bool WC;
};

struct gene
{
  unsigned chrNum;
  string geneName;
  unsigned chrStart;
  unsigned chrEnd;
  bool WC;
  unsigned numTranscripts;
};

struct transcripts{
  unsigned chrNum;
  vector<transcript> transcriptList;
};

struct genes{
  unsigned chrNum;
  vector<gene> geneList;
};

struct exons{
  unsigned chrNum;
  vector<exon> exonList;
};

unsigned get_chr_number(string chrString){
  string::size_type firstGoodPos=chrString.find_first_not_of(" \t\r\n");
  string::size_type lastGoodPos=chrString.find_last_not_of(" \t\r\n");   
  string goodString=chrString.substr(firstGoodPos,lastGoodPos+1-firstGoodPos);
  //cout<<endl<<firstGoodPos;
  //cout<<endl<<lastGoodPos;
  //cout<<endl<<goodString;
  //WAITUSER
  string s_chrNum=goodString.substr(goodString.find_first_not_of("chrCHR"));
  //cout<<endl<<"chr number="<<s_chrNum;
  stringstream ss_chrNum;
  ss_chrNum<<s_chrNum;
  unsigned chrNum;
  ss_chrNum>>chrNum;
  if(chrNum>=1&&chrNum<=22){
    return chrNum;
  }
  else if(chrNum==0){
    if(s_chrNum=="X"){
      return X_CHROMOSOME;
    }
    else if(s_chrNum=="Y"){
      return Y_CHROMOSOME;
    }
    else if(s_chrNum=="M"){
      return M_CHROMOSOME;
    }
  }
  return 0;
}

string get_chr_text(unsigned chrNum){
  if(chrNum>=1&&chrNum<=22){
    stringstream ss_chrNum;
    ss_chrNum<<chrNum;
    string s_chrNum;
    ss_chrNum>>s_chrNum;
    return s_chrNum;
  }
  else if(chrNum==X_CHROMOSOME){
      return "X";
  }
  else if(chrNum==Y_CHROMOSOME){
      return "Y";
  }
  else if(chrNum==M_CHROMOSOME){
      return "M";
  }
  return "NULL";
}

unsigned get_chr_size(unsigned chrNum){
  string chromSizeFileName="chrom.hg19.sizes";
  ifstream chromSizeFile(chromSizeFileName.c_str());
  if(!chromSizeFile.is_open()) {
    cerr<<endl<<"Cannot open"<<chromSizeFileName; WAITUSER;
    return 0;
  }
  string line;
  while(getline(chromSizeFile,line)){
    stringstream ssLine(line);
    string s_chrNum, s_chrSize;
    ssLine>>s_chrNum;
    ssLine>>s_chrSize;
    if(get_chr_number(s_chrNum)==chrNum){
      stringstream ss_chrSize(s_chrSize);
      unsigned i_chrSize=0;
      ss_chrSize>>i_chrSize;
      return i_chrSize;
    }
  }
  return 0;
}

string get_string_from_int(int num){
  ostringstream o_num;
  o_num<<num;
  return o_num.str();
}

struct feature_bitmaps{
  unsigned chrNum;
  vector<bool>featureBitmapChromosome;
};

struct feature_maps{
  unsigned chrNum;
  vector<unsigned short int>featureMapChromosome; //use only if number of data points per chromosome<=65535
  //vector<unsigned>featureMapChromosome; //use only if number of data points per chromosome <=65535
};


unsigned initialize_feature_bitmaps(vector<feature_bitmaps> &featureBitMapsGenome){
  featureBitMapsGenome.resize(M_CHROMOSOME + 1);
  for (unsigned chrNum=1;chrNum<=M_CHROMOSOME;++chrNum){
    if(get_chr_size(chrNum)!=0){
      featureBitMapsGenome[chrNum].featureBitmapChromosome.resize(get_chr_size(chrNum)+1,0);
      //cerr<<endl<<chrNum<<"\t"<<featureBitMapsGenome[chrNum].featureBitmapChromosome.size()-1;
    }
  }
  return featureBitMapsGenome.size();
}

unsigned initialize_feature_maps(vector<feature_maps> &featureMapsGenome){
  featureMapsGenome.resize(M_CHROMOSOME + 1);
  for (unsigned chrNum=1;chrNum<=M_CHROMOSOME;++chrNum){
    if(get_chr_size(chrNum)!=0){
      featureMapsGenome[chrNum].featureMapChromosome.resize(get_chr_size(chrNum)+1,0);
      //cerr<<endl<<chrNum<<"\t"<<featureMapsGenome[chrNum].featureMapChromosome.size()-1;
    }
  }
  return featureMapsGenome.size();
}


ifstream::pos_type get_file_size(string fileName)
{
    ifstream in(fileName.c_str(), ios::ate | ios::binary);
    return in.tellg(); 
}

int main (int argc, char *argv[]){  
  string gencodeFileName="gencode.v19.annotation.gtf";
  ifstream::pos_type fileSize=get_file_size(gencodeFileName);
  if(fileSize<1) {
    cerr<<endl<<"Cannot open"<<gencodeFileName; WAITUSER;
    return 0;
  }
  cerr<<"Gencode annotations file size="<<fileSize/1e9<<" Gigabytes";
  //cerr<<endl<<sizeof(unsigned short int);/*WAITUSER*/
  ifstream gencodeFile;
  gencodeFile.open(gencodeFileName.c_str());
  string line;
  vector<string> lines;
  while(getline(gencodeFile,line)){
    lines.push_back(line);
  }
  vector<genes> genesData;
  genesData.resize(M_CHROMOSOME+1);
  for(unsigned c=1;c<=M_CHROMOSOME;++c){
    genesData[c].chrNum=c;
    struct gene g;
    g.geneName="NULL-GENE";
    genesData[c].geneList.push_back(g);
  }
  vector<transcripts> transcriptsData;
  transcriptsData.resize(M_CHROMOSOME+1);
  for(unsigned c=1;c<=M_CHROMOSOME;++c){
    transcriptsData[c].chrNum=c;
    struct transcript t;
    t.geneName="NULL-TRANSCRIPT";
    transcriptsData[c].transcriptList.push_back(t);
  }
  vector<exons> exonsData;
  exonsData.resize(M_CHROMOSOME+1);
   for(unsigned c=1;c<=M_CHROMOSOME;++c){
    exonsData[c].chrNum=c;
    struct exon e;
    e.geneName="NULL-EXON";
    exonsData[c].exonList.push_back(e);
  }
  
  vector<feature_bitmaps> exonBitmaps;
  vector<feature_maps> promoterMapsW;
  vector<feature_maps> promoterMapsC;
  vector<feature_maps> transcriptMaps;
  initialize_feature_bitmaps(exonBitmaps);
  initialize_feature_maps(transcriptMaps);
  initialize_feature_maps(promoterMapsW);
  initialize_feature_maps(promoterMapsC);
  unsigned nt=0, ng=0, ne=0;
  for(vector<string>::iterator l=lines.begin();l!=lines.end();++l){
      stringstream ssLine(*l);
      vector<string> gtfElements; 
      unsigned countColumns=0;
      while(ssLine.good()){
	string element;
	ssLine>>element;
	++countColumns;
	//cout<<endl<<element<<"\t"<<countColumns;
	gtfElements.push_back(element);
      }
      if(gtfElements.size()>=9){
	string	s_chrNum=gtfElements[0];
	string	gtfSource=gtfElements[1];
	string	gtfFeature=gtfElements[2];
	string	s_chrStart=gtfElements[3];
	string	s_chrEnd=gtfElements[4];
	string	s_gtfScore=gtfElements[5];
	string	gtfStrand=gtfElements[6];
	string	s_gtfFrame=gtfElements[7];
	string	geneName;
	string	geneStatus;
	string 	geneType;
	string	transcriptType;
	string	transcriptName;
	string	transcriptStatus;
	string s_exonNumber;
	for(vector<string>::iterator e=gtfElements.begin() + 8; e!=gtfElements.end(); ++e){
	    if((*e)=="gene_type"){
	      //cout<<endl<<"Gene type "<<*(e+1);//WAITUSER;
	      string s_geneType=*(e+1);
	      string::size_type firstGoodPos=s_geneType.find_first_not_of("\";");
	      string::size_type lastGoodPos=s_geneType.find_last_not_of("\";");   
	      string goodString=s_geneType.substr(firstGoodPos,lastGoodPos+1-firstGoodPos);
	      geneType=goodString;
// 	      cerr<<endl<<geneType;WAITUSER
	      continue;
	    }
	    else if((*e)=="gene_status"){
	      geneStatus=*(e+1);
	      continue;
	    }
    	    else if((*e)=="gene_status"){
	      geneStatus=*(e+1);
	      continue;
	    }
	    else if((*e)=="gene_name"){
	      //cout<<endl<<"Gene name "<<*(e+1);//WAITUSER;
	      string s_geneName=*(e+1);
	      string::size_type firstGoodPos=s_geneName.find_first_not_of("\";");
	      string::size_type lastGoodPos=s_geneName.find_last_not_of("\";");   
	      string goodString=s_geneName.substr(firstGoodPos,lastGoodPos+1-firstGoodPos);
	      geneName=goodString;
	      continue;
	    }
	    else if((*e)=="transcript_type"){
	      //cout<<endl<<"Transcript type "<<*(e+1);//WAITUSER;
	      string s_transcriptType=*(e+1);
	      string::size_type firstGoodPos=s_transcriptType.find_first_not_of("\";");
	      string::size_type lastGoodPos=s_transcriptType.find_last_not_of("\";");   
	      string goodString=s_transcriptType.substr(firstGoodPos,lastGoodPos+1-firstGoodPos);
	      transcriptType=goodString;
	      continue;
	    }
	    else if((*e)=="transcript_status"){
	      //cout<<endl<<"Transcript status "<<*(e+1);
	      transcriptStatus=*(e+1);
	      continue;
	    }
	    else if((*e)=="transcript_name"){
	      //cout<<endl<<"Transcript name "<<*(e+1);//WAITUSER;
	      string s_transcriptName=*(e+1);
	      string::size_type firstGoodPos=s_transcriptName.find_first_not_of("\";");
	      string::size_type lastGoodPos=s_transcriptName.find_last_not_of("\";");   
	      string goodString=s_transcriptName.substr(firstGoodPos,lastGoodPos+1-firstGoodPos);
	      transcriptName=goodString;
	      continue;
	    }
	    else if((*e)=="exon_number"){
	      s_exonNumber=*(e+1);
	    }
	 }
	if(gtfFeature=="gene"&&geneType=="protein_coding"){
	  ++ng;
  	  if((ng%1000)==0){
	    cerr<<endl<<ng<<" genes loaded";
	  }
	  struct gene g;
	  g.chrNum=get_chr_number(s_chrNum);
	  g.geneName=geneName;
	  stringstream ss_chrStart(s_chrStart);
	  ss_chrStart>>g.chrStart;
	  stringstream ss_chrEnd(s_chrEnd);	  
	  ss_chrEnd>>g.chrEnd;
	  g.WC=(gtfStrand=="+"?1:0);
	  genesData[get_chr_number(s_chrNum)].geneList.push_back(g);
	}
	else if(gtfFeature=="transcript"&&geneType=="protein_coding"){//might change later to protein_coding transcript
	 ++nt;
	 struct transcript t;
	  t.chrNum=get_chr_number(s_chrNum);
	  t.geneName=geneName;
	  t.transcriptName=transcriptName;
	  stringstream	ss_chrStart(s_chrStart);
	  ss_chrStart>>t.chrStart;
	  stringstream	ss_chrEnd(s_chrEnd);
	  ss_chrEnd>>t.chrEnd;
	  t.WC=(gtfStrand=="+"?1:0);
	  t.transcriptType=transcriptType;
	  transcriptsData[get_chr_number(s_chrNum)].transcriptList.push_back(t);
	  if(genesData[get_chr_number(s_chrNum)].geneList.back().geneName == geneName){
	    ++genesData[get_chr_number(s_chrNum)].geneList.back().numTranscripts;
	  }
	  for(unsigned chrPos=(t.chrStart); chrPos<=(t.chrEnd); ++chrPos){
	    if(transcriptMaps[t.chrNum].featureMapChromosome[chrPos]==0){
	      transcriptMaps[t.chrNum].featureMapChromosome[chrPos]=(unsigned short int)(transcriptsData[t.chrNum].transcriptList.size()-1);
	    }
	  }
	  if(t.WC==1){
	    if(t.chrStart-PROMOTER_UPSTREAM <1 || t.chrStart+PROMOTER_DOWNSTREAM >= promoterMapsW[t.chrNum].featureMapChromosome.size()){
	      cerr<<endl<<"Problem allocating promoter region of gene"<<t.geneName<<"\t"<<t.chrNum<<"\t"<<promoterMapsW[t.chrNum].featureMapChromosome.size()<<"\t"<<t.chrStart;//WAITUSER
	    }
	    else
	    for(unsigned chrPos=(t.chrStart)-PROMOTER_UPSTREAM; chrPos<=(t.chrStart)+PROMOTER_DOWNSTREAM; ++chrPos){
	      if(promoterMapsW[t.chrNum].featureMapChromosome[chrPos]==0){
		promoterMapsW[t.chrNum].featureMapChromosome[chrPos]=(unsigned short int)(transcriptsData[t.chrNum].transcriptList.size()-1);
	      }
	    }
	  }
	  else{
	    if(t.chrEnd+PROMOTER_UPSTREAM >=promoterMapsC[t.chrNum].featureMapChromosome.size() || t.chrEnd-PROMOTER_DOWNSTREAM<1){
	      cerr<<endl<<"Problem allocating promoter region of gene"<<t.geneName<<"\t"<<t.chrNum<<"\t"<<promoterMapsC[t.chrNum].featureMapChromosome.size()<<"\t"<<t.chrStart;//WAITUSER
	    }
	    else
	    for(unsigned chrPos=(t.chrEnd)+PROMOTER_UPSTREAM; chrPos>=(t.chrEnd)-PROMOTER_DOWNSTREAM; --chrPos){
	      if(promoterMapsC[t.chrNum].featureMapChromosome[chrPos]==0){
		promoterMapsC[t.chrNum].featureMapChromosome[chrPos]=(unsigned short int)(transcriptsData[t.chrNum].transcriptList.size()-1);
	      }
	    }
	  }
	      
	}
	else if(gtfFeature=="exon"&&geneType=="protein_coding"){
	  ++ne;
	  struct exon e;
	  e.chrNum=get_chr_number(s_chrNum);
	  e.geneName=geneName;
	  stringstream ss_exonNumber(s_exonNumber);
	  ss_exonNumber>>e.exonNumber;
	  e.transcriptName=transcriptName;
	  stringstream ss_chrStart(s_chrStart);
	  ss_chrStart>>e.chrStart;
	  stringstream ss_chrEnd(s_chrEnd);
	  ss_chrEnd>>e.chrEnd;
	  e.WC=(gtfStrand=="+"?1:0);
	  //cerr<<endl<<"Exon number:"<<s_exonNumber;
	  exonsData[get_chr_number(s_chrNum)].exonList.push_back(e);
	  for(unsigned chrPos=(e.chrStart); chrPos<=(e.chrEnd); ++chrPos){
	    exonBitmaps[e.chrNum].featureBitmapChromosome[chrPos]=1;
	  }
	}

      }
    //WAITUSER
  }
cerr<<endl<<lines.size();
cerr<<endl<<"num transcripts="<<nt<<"\t num genes="<<ng<<"\t num exons="<<ne;
for (unsigned c=1;c<=M_CHROMOSOME;++c){
  if(transcriptsData[c].transcriptList.size()>1){
    cerr<<endl<<c<<"\t"<<transcriptsData[c].transcriptList.size()<<"\t"<<genesData[c].geneList.size();
  }
}
lines.clear();
{
  string peakFileName="";
  if(argc>1){
    peakFileName=argv[1];
  }
  else{
    cerr<<endl<<"Input should be a 10-column file with ERG peaks (summit location) in the 1st 5 columns and"<< \
		" co-occupying EZH2 peaks (summit location) in the next 5 columns"<< endl << \
		"chr1    762714  762715  MACS_peak_7     439.21  chr1    763119  763120  MACS_peak_1     58.87";
    peakFileName="test.bed";
  }
  ifstream::pos_type fileSize=get_file_size(peakFileName);
  cerr<<endl<<"Using input file"<<peakFileName;
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
      //cerr<<endl<<"Peak Mid="<<peakMid;
      string s_chrNum2=".";
      string s_chrStart2=".";
      string s_chrEnd2=".";
      unsigned chrNum2=0;
      unsigned chrStart2=0;
      unsigned chrEnd2=0;
      unsigned peakMid2=0;
      string name2=".";
      string s_score2=".";
      if(bedElements.size()==10){//each line has EZH2 data as well
	s_chrNum2=bedElements[5];
	s_chrStart2=bedElements[6];
	s_chrEnd2=bedElements[7];
	name2=bedElements[8];
	s_score2=bedElements[9];
	chrNum2=get_chr_number(s_chrNum2);
	stringstream ss_chrStart2(s_chrStart2);
	ss_chrStart2>>chrStart2;
	stringstream ss_chrEnd2(s_chrEnd2);
	ss_chrEnd2>>chrEnd2;
	if(chrNum2!=chrNum){
	  cerr<<endl<<"ERROR in BED Line"<<*l;WAITUSER;
	}
	peakMid2=chrStart2/2+chrEnd2/2;
      }
      int transcriptStartSite=-1;
      unsigned geneWC=0;
      unsigned short int promoterIndexTestW=promoterMapsW[chrNum].featureMapChromosome[peakMid];
      unsigned short int promoterIndexTestC=promoterMapsC[chrNum].featureMapChromosome[peakMid];
      unsigned short int transcriptIndexTest=transcriptMaps[chrNum].featureMapChromosome[peakMid];
      if(promoterIndexTestW>0||promoterIndexTestC>0){
	if(promoterIndexTestW>0){
	  cout<<endl<<"chr"<<get_chr_text(chrNum)<<"\t"<<peakMid<<"\t"<<name<<"\t"<<score;
	  cout<<"\tPROMOTER";
	  cout<<"\t"<<transcriptsData[chrNum].transcriptList[promoterIndexTestW].geneName;
	  geneWC=transcriptsData[chrNum].transcriptList[promoterIndexTestW].WC;
	  cout<<"\t"<<(geneWC==1?"+":"-");
	  cout<<"\t"<<transcriptsData[chrNum].transcriptList[promoterIndexTestW].transcriptName;
	  cout<<"\t"<<transcriptsData[chrNum].transcriptList[promoterIndexTestW].transcriptType;
	  transcriptStartSite=(transcriptsData[chrNum].transcriptList[promoterIndexTestW].WC==1? \
			      transcriptsData[chrNum].transcriptList[promoterIndexTestW].chrStart: \
			      transcriptsData[chrNum].transcriptList[promoterIndexTestW].chrEnd); 
	  cout<<"\t"<<transcriptStartSite;
	  cout<<"\t"<<(transcriptsData[chrNum].transcriptList[promoterIndexTestW].WC==1? \
		      (int)peakMid-(int)transcriptsData[chrNum].transcriptList[promoterIndexTestW].chrStart: \
		      (int)transcriptsData[chrNum].transcriptList[promoterIndexTestW].chrEnd-(int)peakMid);
	    if(transcriptStartSite > -1 ){
	      cout<<"\t"<<peakMid2<<"\t"<<name2<<"\t"<<s_score2<<"\t"<<(int)(peakMid2-peakMid) \
	      <<"\t"<<(geneWC==1?(int)peakMid2-(int)transcriptStartSite:(int)transcriptStartSite-(int)peakMid2);
	    }
	    else{
	      cout<<"\t"<<peakMid2<<"\t"<<name2<<"\t"<<s_score2<<"\t"<<(int)(peakMid2-peakMid)<<"\t.";
	    }
	}
	if(promoterIndexTestC>0){
	  cout<<endl<<"chr"<<get_chr_text(chrNum)<<"\t"<<peakMid<<"\t"<<name<<"\t"<<score;
	  cout<<"\tPROMOTER";
	  cout<<"\t"<<transcriptsData[chrNum].transcriptList[promoterIndexTestC].geneName;
	  geneWC=transcriptsData[chrNum].transcriptList[promoterIndexTestC].WC;
	  cout<<"\t"<<(geneWC==1?"+":"-");
	  cout<<"\t"<<transcriptsData[chrNum].transcriptList[promoterIndexTestC].transcriptName;
	  cout<<"\t"<<transcriptsData[chrNum].transcriptList[promoterIndexTestC].transcriptType;
	  transcriptStartSite=(transcriptsData[chrNum].transcriptList[promoterIndexTestC].WC==1? \
			      transcriptsData[chrNum].transcriptList[promoterIndexTestC].chrStart: \
			      transcriptsData[chrNum].transcriptList[promoterIndexTestC].chrEnd); 
	  cout<<"\t"<<transcriptStartSite;
	  cout<<"\t"<<(transcriptsData[chrNum].transcriptList[promoterIndexTestC].WC==1? \
		      (int)peakMid-(int)transcriptsData[chrNum].transcriptList[promoterIndexTestC].chrStart: \
		      (int)transcriptsData[chrNum].transcriptList[promoterIndexTestC].chrEnd-(int)peakMid);
  	  if(transcriptStartSite > -1 ){
	    cout<<"\t"<<peakMid2<<"\t"<<name2<<"\t"<<s_score2<<"\t"<<(int)(peakMid2-peakMid) \
	    <<"\t"<<(geneWC==1?(int)peakMid2-(int)transcriptStartSite:(int)transcriptStartSite-(int)peakMid2);
	  }
	  else{
	    cout<<"\t"<<peakMid2<<"\t"<<name2<<"\t"<<s_score2<<"\t"<<(int)(peakMid2-peakMid)<<"\t.";
	  }
	}
      }
      else if(transcriptIndexTest>0){
	  cout<<endl<<"chr"<<get_chr_text(chrNum)<<"\t"<<peakMid<<"\t"<<name<<"\t"<<score;
	  bool exonTest=0;
	  exonTest=exonBitmaps[chrNum].featureBitmapChromosome[peakMid];
	  if(exonTest){
	    cout<<"\tEXON";
	  }
	  else{
	    cout<<"\tINTRON";
	  }
	  cout<<"\t"<<transcriptsData[chrNum].transcriptList[transcriptIndexTest].geneName;
	  geneWC=transcriptsData[chrNum].transcriptList[transcriptIndexTest].WC;
	  cout<<"\t"<<(geneWC==1?"+":"-");
	  cout<<"\t"<<transcriptsData[chrNum].transcriptList[transcriptIndexTest].transcriptName;
  	  cout<<"\t"<<transcriptsData[chrNum].transcriptList[transcriptIndexTest].transcriptType;
	  transcriptStartSite=(transcriptsData[chrNum].transcriptList[transcriptIndexTest].WC==1? \
			      transcriptsData[chrNum].transcriptList[transcriptIndexTest].chrStart: \
			      transcriptsData[chrNum].transcriptList[transcriptIndexTest].chrEnd);	
	  cout<<"\t"<<transcriptStartSite;
	  cout<<"\t"<<(transcriptsData[chrNum].transcriptList[transcriptIndexTest].WC==1? \
		      (int)peakMid-(int)transcriptsData[chrNum].transcriptList[transcriptIndexTest].chrStart: \
		      (int)transcriptsData[chrNum].transcriptList[transcriptIndexTest].chrEnd-(int)peakMid);
  	  if(transcriptStartSite > -1 ){
	    cout<<"\t"<<peakMid2<<"\t"<<name2<<"\t"<<s_score2<<"\t"<<(int)(peakMid2-peakMid) \
	    <<"\t"<<(geneWC==1?(int)peakMid2-(int)transcriptStartSite:(int)transcriptStartSite-(int)peakMid2);
	  }
	  else{
	    cout<<"\t"<<peakMid2<<"\t"<<name2<<"\t"<<s_score2<<"\t"<<(int)(peakMid2-peakMid)<<"\t.";
	  }
      }
      else{
	cout<<endl<<"chr"<<get_chr_text(chrNum)<<"\t"<<peakMid<<"\t"<<name<<"\t"<<score;
	cout<<"\tINTERGENIC\t.\t.\t.\t.\t.\t.";
	if(transcriptStartSite > -1 ){
	    cout<<"\t"<<peakMid2<<"\t"<<name2<<"\t"<<s_score2<<"\t"<<(int)(peakMid2-peakMid) \
	    <<"\t"<<(geneWC==1?(int)peakMid2-(int)transcriptStartSite:(int)transcriptStartSite-(int)peakMid2);
	}
	else{
	    cout<<"\t"<<peakMid2<<"\t"<<name2<<"\t"<<s_score2<<"\t"<<(int)(peakMid2-peakMid)<<"\t.";
	}
      }

    }
  }
  lines.clear();
}
  exonBitmaps.clear();
  promoterMapsW.clear();
  promoterMapsC.clear();
  transcriptMaps.clear();
  return 0;
  
}
