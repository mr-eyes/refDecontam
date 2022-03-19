#include "kDataFrame.hpp"
#include <stdexcept>
#include "tuple"
#include <sys/stat.h>
#include "colored_kDataFrame.hpp"
#include <kseq/kseq.h>
#include <zlib.h>
#include <cstdio>
#include <iostream>
#include <vector>
#include <glob.h> // glob(), globfree()
#include <cassert>
#include "kmerDecoder.hpp"
#include <string>     // std::string, std::to_string

#define KSIZE 21
#define CHUNK_SIZE 5000

using namespace std;
using namespace phmap;


string create_dir(string output_file, int serial) {
    int dir_err;
    string new_name = "";

    if (!serial) {
        dir_err = mkdir(output_file.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        new_name = output_file;
    }
    else {
        new_name = output_file + "_v." + std::to_string(serial);
        dir_err = mkdir(new_name.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
    }

    if (-1 == dir_err) return create_dir(output_file, ++serial);

    return new_name;
}

class fileHandler {

public:
    ofstream fileStream;

    fileHandler(string& filename) {
        this->fileStream.open(filename);
    }

    void write(string& line) {
        this->fileStream << line;
    }

    void close() {
        fileStream.close();
    }

};

inline bool file_exists(const std::string& name) {
    struct stat buffer
    {
    };
    return (stat(name.c_str(), &buffer) == 0);
}

inline string time_diff(std::chrono::high_resolution_clock::time_point& t1) {
    std::chrono::high_resolution_clock::time_point t2 = std::chrono::high_resolution_clock::now();
    auto milli = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
    long hr = milli / 3600000;
    milli = milli - 3600000 * hr;
    long min = milli / 60000;
    milli = milli - 60000 * min;
    long sec = milli / 1000;
    milli = milli - 1000 * sec;
    string timeDiff;
    timeDiff.append(to_string(min));
    timeDiff.append(":");
    timeDiff.append(to_string(sec));
    timeDiff.append(":");
    timeDiff.append(to_string(milli));

    return timeDiff;
}

// thanks to https://stackoverflow.com/a/8615450/3371177
std::vector<std::string> glob(const std::string& pattern) {
    using namespace std;

    // glob struct resides on the stack
    glob_t glob_result;
    memset(&glob_result, 0, sizeof(glob_result));

    // do the glob operation
    int return_value = glob(pattern.c_str(), GLOB_TILDE, NULL, &glob_result);
    if (return_value != 0) {
        globfree(&glob_result);
        stringstream ss;
        ss << "glob() failed with return_value " << return_value << endl;
        throw std::runtime_error(ss.str());
    }

    // collect all the filenames into a std::list<std::string>
    vector<string> filenames;
    for (size_t i = 0; i < glob_result.gl_pathc; ++i) {
        filenames.push_back(string(glob_result.gl_pathv[i]));
    }

    // cleanup
    globfree(&glob_result);

    // done
    return filenames;
}

tuple<string, vector<int>> score(vector<uint32_t>& genomes) {

    vector<int> sources;

    if (genomes.empty())
        return make_tuple("unmapped", sources);

    flat_hash_map<int, int> scores;
    flat_hash_map<int, int> reverse_scores;
    flat_hash_map<int, int> countFreq;
    vector<int> all_scores;

    for (const auto& genome : genomes) {
        scores[genome]++;
    }

    for (const auto& score : scores) {
        countFreq[score.second]++;
        all_scores.emplace_back(score.second);
        reverse_scores[score.second] = score.first;
    }

    auto max = std::max_element(all_scores.begin(), all_scores.end());

    if (countFreq[*max] == 1) {
        sources.emplace_back(reverse_scores[*max]);
        return make_tuple("unique", sources);
    }

    for (const auto& score : scores) {
        if (score.second == *max) {
            sources.emplace_back(score.first);
        }
    }

    return make_tuple("ambig", sources);

}



int main(int argc, char** argv) {

    KS_FULL_COMMENT = true;

    if (argc != 3) {
        cerr << "run ./decontaminate <ref_fasta> <index_prefix>" << endl;
        exit(1);
    }

    string reads_file = argv[1];
    string kfs_dir = argv[2];

    kDataFrame* frame;
    std::string dir_prefix = kfs_dir.substr(kfs_dir.find_last_of("/\\") + 1);

    flat_hash_map<string, string> namesMap;
    string names_fileName = kfs_dir;

    flat_hash_map<string, uint64_t> tagsMap;
    flat_hash_map<string, uint64_t> groupNameMap;
    flat_hash_map<uint64_t, std::vector<uint32_t>>* legend = new flat_hash_map<uint64_t, std::vector<uint32_t>>();
    flat_hash_map<uint64_t, uint64_t> colorsCount;
    uint64_t readID = 0, groupID = 1;
    string seqName, groupName;
    string line;

    int total_kfs_number = 0;

    // get kSize and type
    for (const auto& dirEntry : glob(kfs_dir + "/*")) {
        string file_name = (string)dirEntry;
        size_t lastindex = file_name.find_last_of(".");
        string kf_prefix = file_name.substr(0, lastindex);
        std::string::size_type idx;
        idx = file_name.rfind('.');
        std::string extension = "";
        if (idx != std::string::npos) extension = file_name.substr(idx + 1);
        int detected_kSize;
        hashingModes _hm;
        if (extension == "mqf" || extension == "phmap") {
            auto* _kf = kDataFrame::load(kf_prefix);
            _hm = _kf->getkmerDecoder()->hash_mode;
            detected_kSize = _kf->getkSize();
            cout << "Detected kSize: " << detected_kSize << endl;
        }
        else {
            continue;
        }
        // if(extension == "mqf") {frame = new kDataFrameMQF(detected_kSize, 30, mumur_hasher); break;} // temp. switch off
        if (extension == "mqf") { frame = new kDataFramePHMAP(detected_kSize, _hm); break; }
        else if (extension == "phmap") { frame = new kDataFramePHMAP(detected_kSize, _hm); break; }
        else { continue; }
    }

    cout << "namesmap construction done..." << endl;

    if (!file_exists(reads_file)) {
        throw std::runtime_error("Could not open the reads fasta file");
    }

    flat_hash_map<string, kDataFrame*> all_kfs;
    flat_hash_map<string, int> genome_to_id;
    flat_hash_map<int, string> id_to_genome;



    // Loading kfs
    int counter = 0;
    for (const auto& dirEntry : glob(kfs_dir + "/*")) {
        string file_name = (string)dirEntry;
        size_t lastindex = file_name.find_last_of(".");
        string kf_prefix = file_name.substr(0, lastindex);

        std::string kf_basename = kf_prefix.substr(kf_prefix.find_last_of("/\\") + 1);


        std::string::size_type idx;
        idx = file_name.rfind('.');
        std::string extension = "";
        if (idx != std::string::npos) extension = file_name.substr(idx + 1);
        if (extension != "mqf" and extension != "phmap") continue;

        all_kfs[kf_basename] = kDataFrame::load(kf_prefix);
        genome_to_id[kf_basename] = counter++;
        id_to_genome[genome_to_id[kf_basename]] = kf_basename;
    }

    // kProcessor Index Loading
    map<string, fileHandler*> fasta_writer;

    cerr << "Creating fasta file handlers..." << endl;
    for (auto& item : all_kfs) {
        string file_name = "genome_" + item.first + "_partition.fa";
        fasta_writer[item.first] = new fileHandler(file_name);
    }

    string unmapped_file_name = "unmapped_partition.fa";
    fasta_writer["unmapped"] = new fileHandler(unmapped_file_name);

    int chunkSize = CHUNK_SIZE;
    int kSize = KSIZE;


    gzFile fp;
    kseq_t* kseqObj;
    fp = gzopen(reads_file.c_str(), "r");
    kseqObj = kseq_init(fp);

    cout << "Processing started ..." << endl;

    auto* KD = kmerDecoder::getInstance(KMERS, mumur_hasher, { {"kSize", kSize} });


    int total = 0;
    int chunks = 0;

    uint64_t _unique = 0;
    uint64_t _ambig = 0;
    uint64_t _unmatched = 0;

    while (kseq_read(kseqObj) >= 0) {
        // if (string(kseqObj->seq.s).size() < kSize) continue;

        std::string seq = kseqObj->seq.s;
        std::string id;
        id.append(kseqObj->name.s);
        if (kseqObj->comment.l) id.append(kseqObj->comment.s);

        std::string record = ">";
        record.append(id);
        record.append("\n");
        record.append(seq);
        record.append("\n");

        vector<uint32_t> kmers_matches;

        for (auto const& ref : all_kfs) {
            for (unsigned long i = 0; i < seq.size() - kSize + 1; i++) {
                uint64_t kmer = KD->hash_kmer(seq.substr(i, kSize));
                if (ref.second->getCount(kmer)) {
                    kmers_matches.emplace_back(genome_to_id[ref.first]);
                }
            }
        }

        auto category = score(kmers_matches);


        if (get<0>(category) == "unique") {
            _unique++;
            fasta_writer[id_to_genome[get<1>(category)[0]]]->write(record);
        }
        else if (get<0>(category) == "unmapped") {
            _unmatched++;
            fasta_writer["unmapped"]->write(record);
        }
        else if (get<0>(category) == "ambig") {
            _ambig++;
            for (auto const& genomeID : get<1>(category)) {
                fasta_writer[id_to_genome[genomeID]]->write(record);
            }
        }


        total++;
        if (total == 5000) {
            cout << "processed: " << 5000 * ++chunks << "contigs" << endl;
            total = 0;
        }

    }

    cout << "_unique: " << _unique << endl;
    cout << "_ambig: " << _ambig << endl;
    cout << "_unmatched: " << _unmatched << endl;

    for (auto f : fasta_writer)
        f.second->close();


    return 0;
}