#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <math.h>
#include <cstring>
using namespace std;


const int MAX_FASTA_READING = 50;

int fastawrite;
ofstream logFile;


class Data {
public:
	string chrom;
	int pos;
	string ref;
	string alt;


	Data() {
		this->chrom = "";
		this->pos = -1;
		this->ref = "";
		this->alt = "";
	}
	void print() { logFile << this->chrom << "\t" << this->pos << "\t" << this->ref << "\t" << this->alt << endl; }
	bool validate()
	{
		if (this->chrom == "" || this->pos == -1 || this->ref == "" || this->alt == "") return false;
		else if ( this->ref.length() > 1 || this->alt.length() > 1 )	return false;
		else return true;
	}
};
void VCF_read(ifstream & input, vector<Data>& vcf, vector<string>& scaffolds)
{
	while (!input.eof())
	{
		string tmp_VCF;
		do { getline(input, tmp_VCF);

			if (tmp_VCF[0] == '#' && tmp_VCF[2] == 'c')
			{
				int tPos = tmp_VCF.find_first_of('=', 0);
				string scaffold = tmp_VCF.substr(tPos + 1);

				int sPos = scaffold.find_first_of('=', 0)+1;
				int ePos = scaffold.find_first_of(',', 0)-sPos;

				scaffold = scaffold.substr(sPos, ePos);
				scaffolds.push_back(scaffold);
			}

		} while (tmp_VCF[0] == '#');


		istringstream vcfStream(tmp_VCF);
		string token;
		int format = 0;

		Data change;
		while (vcfStream >> token)
		{
			if (format == 0)
			{
				change.chrom = token;
			}
			else if (format == 1)
            {
                change.pos = atoi( token.c_str() );
                // change.pos = stoi(token);
            }
			else if (format == 3) change.ref = token;
			else if (format == 4) change.alt = token;
			else {}

			format++;
		}

		if (change.validate() == true) vcf.push_back(change);
		else {}
	}
}
int fa_Read(ifstream & input, string & seq)
{
	if (input.eof()) return -1;

	getline(input, seq);
	if (seq[0] == '>') getline(input, seq);
	return seq.length();
}
void fa_Write(ofstream & output, string seq)
{
	output << seq << endl;
}
int CheckVCFRange(const vector<Data> input, int now, vector<Data> & output)
{
	int pass = now + 1;
	logFile << "Next vcf Index ... \t" << pass << endl;
	string chrom = output.back().chrom;
	int pos = output.back().pos;


	int line = (pos -1) / MAX_FASTA_READING;
	logFile << "Now vcf line ... \t" << line << endl;

	for (pass; pass < input.size(); pass++)
	{
		if (chrom != input.at(pass).chrom) return pass;
		if ((input.at(pass).pos -1)% MAX_FASTA_READING == 0) return pass;

		if (line == ((input.at(pass).pos -1) / MAX_FASTA_READING)) {
			logFile << "input vcf Pos : " << input.at(pass).pos << endl;
			output.push_back(input.at(pass));
		}
		else {
			logFile << "line Last Inputed vcf Pos = " << output.back().pos << endl;
			logFile << "Output vcf data index ... \t" << pass << endl;
			return pass;
		}
	}
	return pass;
}
int seq_check(const vector<Data> vcfRange)
{
	Data change = vcfRange.back();
	double seqLen = (double)MAX_FASTA_READING;

	int vPos = change.pos % (int)seqLen;
	int vLength = change.ref.length();
	double refRange = vPos + vLength;
	seqLen = ceil(refRange / seqLen);

	return seqLen - 1;
}
void fa_change(vector<Data> vcfRange, string & seq)
{
	logFile << "change count " << vcfRange.size() << endl;
	int vPos;
	int preLen = 0;
	int split = 0;

	string changeSeq = "";

	for (int c=0; c< vcfRange.size(); c++)
    {
        Data change = vcfRange.at(c);
		vPos = change.pos % MAX_FASTA_READING;
		preLen = vPos - 1;
		if (vPos == 0)  preLen = 49;

		change.alt = change.alt[0];

		string pre = seq.substr(split, preLen - split);
		logFile << "pre Seq : " << pre << " + "
			<< "change Seq : " << change.ref << " -> " << change.alt << endl;
		changeSeq += (pre + change.alt);
		logFile << "changed Seq  : " << changeSeq << endl;
		//changeSeq = pre + alt + pro
		split = preLen + change.ref.length();
    }
/*
	for each ( Data change in vcfRange)
	{
		logFile << "change Pos VCF : " << change.pos << endl;
		vPos = change.pos % MAX_FASTA_READING;
		preLen = vPos - 1;
		if (vPos == 0)  preLen = 49;

		change.alt = change.alt[0];

		string pre = seq.substr(split, preLen - split);
		logFile << "pre Seq : " << pre << " + "
			<< "change Seq : " << change.ref << " -> " << change.alt << endl;
		changeSeq += (pre + change.alt);
		logFile << "changed Seq  : " << changeSeq << endl;
		//changeSeq = pre + alt + pro
		split = preLen + change.ref.length();
	}
*/
	changeSeq += seq.substr(split);
	seq = changeSeq;
	vcfRange.clear();
}


int main(int argc, char* argv[])
{
    if( argc < 4 )
    {
        cout << " Parameter 1 - VCF file ( .vcf file )" << endl;
        cout << " Parameter 2 - Fasta files directory ( Must be split by chromosome name )" << endl;
        cout << " Parameter 3 - Output prefix" << endl;

        return 0;
    }

    char* vcf_Dir = argv[1]; //  char* vcf_Dir = "E:/YEAST/vcf/sk1.masurca.reverse/sk1.masurca.snps.vcf";
    char* fasta_Dir = argv[2];   //  char* fasta_Dir = "E:/YEAST/vcf/sk1.masurca.reverse/byName/";
    char* prefix = argv[3];

      char* out_with_extension = (char*)malloc( strlen(prefix) + 1 + 8 );
      char* log_with_extension = (char*)malloc( strlen(prefix) + 1 + 4 );

    strcpy(out_with_extension, prefix);
    char* out_Dir = strcat(out_with_extension, ".edit.fa"); //  char* out_Dir = "E:/YEAST/vcf/sk1.masurca.reverse.edit.fa";
    strcpy(log_with_extension, prefix);
    char* log_DIr = strcat(log_with_extension, ".log");    //   char* log_DIr = "E:/YEAST/vcf/sk1.masurca.reverse/sk1.masurca.reverse.log.txt";



	ifstream vcf_File;
	ifstream fa_File;
	ofstream fa_OFile;

	vector<Data> vcf;
	vector<string> scaffolds;

	cout << "start" << endl;

	vcf_File.open(vcf_Dir);
	fa_OFile.open(out_Dir);
	logFile.open(log_DIr);

	logFile << "vcf reading ...";
	VCF_read(vcf_File, vcf, scaffolds);
	logFile << "end !" << endl;

	string CHROM_FA = "Init";
	string seq;
	int faLine;

	cout << scaffolds.size() << endl;

	for (int i = 0; i < vcf.size(); )
	{
		vector <Data> line;

		if (CHROM_FA != vcf.at(i).chrom)
		{
			cout << CHROM_FA << ",\t" << scaffolds.front() << endl;

			logFile << CHROM_FA << "\t finished...";
			if (fa_File.is_open()) { while (fa_Read(fa_File, seq) != -1) fa_Write(fa_OFile, seq); }
			fa_File.clear();
			fa_File.close();
			logFile << "\t Write rest of " << CHROM_FA << endl;
			CHROM_FA = vcf.at(i).chrom;

			while (CHROM_FA != scaffolds.front())
			{
				fa_OFile << ">" << scaffolds.front() << endl;
				string faName = fasta_Dir + scaffolds.front() + ".fa";
				fa_File.open(faName.c_str());
				while (fa_Read(fa_File, seq) != -1) fa_Write(fa_OFile, seq);

				scaffolds.erase(scaffolds.begin());
				fa_File.close();
				cout << CHROM_FA << ",\t" << scaffolds.front() << endl;

			}
			fa_OFile << ">" << CHROM_FA << endl;
			logFile << "Open " << CHROM_FA << "..." << endl;
            string faName = fasta_Dir + CHROM_FA + ".fa";
			fa_File.open(faName.c_str());
			fastawrite = 0;
			faLine = 0;
			scaffolds.erase( remove( scaffolds.begin(), scaffolds.end(), vcf.at(i).chrom ), scaffolds.end() );
		}
		line.push_back(vcf.at(i));
		logFile << "pre -check line_push_back Data pos : " << vcf.at(i).pos << endl;
		int vLine = (vcf.at(i).pos - 1) / MAX_FASTA_READING;

		while (faLine < vLine)
		{
			fa_Read(fa_File, seq);	faLine++;
			fa_Write(fa_OFile, seq);
		}
		fa_Read(fa_File, seq); faLine++;

		while (true)
		{
			//save 0~80 range changes;
			logFile << "Check VCF Range " << endl;
			i = CheckVCFRange(vcf, i, line);
			logFile << "index " << i << " :: " << vcf.size() << endl;
			if (i == vcf.size()) break;
			logFile << "input the change data to line vector ..." << endl;
			// seq check may be refenence can be up to Length MAX_FASTA_READING

			int limit = seq_check(line);
			logFile << "Extend sequence : \t" << limit << endl;
			for (int j = 0; j < limit; j++)
			{
				string more = "";
				fa_Read(fa_File, more);	faLine++;
				seq += more;
			}

			logFile << (vcf.at(i).pos -1) / MAX_FASTA_READING << "\t:\t" << (line.back().pos -1) / MAX_FASTA_READING + limit << endl;
			if (((vcf.at(i).pos -1) / MAX_FASTA_READING) <= (((line.back().pos -1) / MAX_FASTA_READING) + limit)) {

				if (vcf.at(i).chrom != line.back().chrom) {
					logFile << vcf.at(i).chrom << "\t:\t" << line.back().chrom << endl;
					break;
				}
				line.push_back(vcf.at(i));
				logFile << "pre -check line_push_back Data pos : " << vcf.at(i).pos << endl;

			}
			else {
				logFile << "Decide the sequence Length to change ... \t " << seq.length() << endl;
				break;
			}
		}

		logFile << "fasta change " << endl;
		logFile << "Unchanged Sequence :\t" << seq << endl;
		fa_change(line, seq);
		logFile << "changed Sequence :\t" << seq << endl;

		logFile << "fa_Write ..." << endl;
		fa_Write(fa_OFile, seq);
		logFile << endl;
	}

	if (fa_File.is_open()) { while (fa_Read(fa_File, seq) != -1) fa_Write(fa_OFile, seq); }

	for (int i = 0; i < scaffolds.size(); i++)
	{
		fa_File.close();
		string faName = fasta_Dir + scaffolds.at(i);
		fa_OFile << '>' << scaffolds.at(i) << endl;
		faName += ".fa";
		fa_File.open(faName.c_str());

		while (fa_Read(fa_File, seq) != -1) fa_Write(fa_OFile, seq);
		logFile << "write " << faName << " scaffold " << i << endl;
	}

	fa_File.close();
	fa_OFile.close(); logFile.close();

	return 0;
}

